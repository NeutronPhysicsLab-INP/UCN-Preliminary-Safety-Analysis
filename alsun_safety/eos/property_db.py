"""
NIST-tabulated equation-of-state handler for liquid / gaseous deuterium.

Loads the saturation line and a stack of isothermal tables, then provides
:meth:`PropertyDB.state_from_T_rho` to evaluate the fluid state (P, u, h,
cp, cv, k, gamma, phase) at any (T, rho) inside the loaded grid.

The phase-detection logic implements the algorithm in Section 3.1 of the
manuscript:

    1. Compare specific volume v = 1/rho to the saturated bounds v_l(T),
       v_g(T) from the saturation table.
    2. If v_l < v < v_g, treat as two-phase: clamp P to P_sat(T), use the
       lever rule for extensive properties.
    3. Otherwise treat as single-phase and bilinearly interpolate on the
       isothermal tables (rho -> P -> {u,h,cp,cv,k}).

Original implementation by the AlSUN team; refactored into a standalone
module for the public release. See ``alsun_safety/eos/solver.py`` for
inversion utilities (T from u, rho from P).
"""

from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Optional, Tuple
import os
import re

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Numeric utilities
# ---------------------------------------------------------------------------

def _clamp(x: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, x))


def _interp1(xq: np.ndarray, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """1-D linear interpolation with clamping at the ends."""
    x  = np.asarray(x,  float)
    y  = np.asarray(y,  float)
    xq = np.asarray(xq, float)
    xq_clamped = np.clip(xq, x[0], x[-1])
    return np.interp(xq_clamped, x, y)


def _unique_monotone(x: np.ndarray, y: np.ndarray, increasing: bool = True
                     ) -> Tuple[np.ndarray, np.ndarray]:
    """Sort by x, drop duplicates, return strictly monotone x with aligned y."""
    x = np.asarray(x, dtype=float).copy()
    y = np.asarray(y, dtype=float).copy()
    idx = np.argsort(x)
    x, y = x[idx], y[idx]
    x_u, idx_u = np.unique(x, return_index=True)
    x, y = x[idx_u], y[idx_u]
    if increasing:
        mask = np.concatenate([[True], np.diff(x) > 0])
        return x[mask], y[mask]
    # decreasing case
    x, y = -x, y
    idx = np.argsort(x)
    x, y = x[idx], y[idx]
    x_u, idx_u = np.unique(x, return_index=True)
    x, y = x[idx_u], y[idx_u]
    mask = np.concatenate([[True], np.diff(x) > 0])
    x, y = x[mask], y[mask]
    return -x, y


def _to_num(series: pd.Series, factor: float = 1.0) -> np.ndarray:
    """Coerce NIST table strings to numeric; map 'undefined', dashes, etc. to NaN."""
    s = series.astype(str).str.strip().replace({
        "undefined": np.nan, "UNDEFINED": np.nan,
        "nan": np.nan, "NaN": np.nan, "N/A": np.nan, "NA": np.nan,
        "—": np.nan, "-": np.nan, "": np.nan,
    })
    return (pd.to_numeric(s, errors="coerce") * factor).to_numpy()


def _interp_fill_by_T(T: np.ndarray, arr: np.ndarray) -> np.ndarray:
    """Fill NaNs along T using linear interpolation; clamp at ends."""
    T = np.asarray(T, float)
    y = np.asarray(arr, float)
    mask = np.isfinite(y)
    if mask.sum() == 0:
        return np.full_like(y, np.nan)
    if mask.sum() == 1:
        return np.full_like(y, y[mask][0])
    y2 = y.copy()
    y2[~mask] = np.interp(T[~mask], T[mask], y[mask])
    return y2


# ---------------------------------------------------------------------------
# Saturation properties
# ---------------------------------------------------------------------------

@dataclass
class SatProps:
    """Tabulated saturation-line properties for a single fluid."""
    T_grid: np.ndarray
    P_sat:  np.ndarray
    v_l:    np.ndarray
    v_g:    np.ndarray
    u_l:    np.ndarray
    u_g:    np.ndarray
    h_l:    np.ndarray
    h_g:    np.ndarray
    cp_l:   np.ndarray
    cp_g:   np.ndarray
    cv_l:   np.ndarray
    cv_g:   np.ndarray
    k_l:    np.ndarray
    k_g:    np.ndarray

    def eval(self, T: float) -> Dict[str, float]:
        """Interpolate all saturation properties at temperature T."""
        Tq = float(_clamp(T, float(self.T_grid.min()), float(self.T_grid.max())))
        it = lambda arr: float(_interp1(np.array([Tq]), self.T_grid, arr)[0])
        return dict(
            P_sat=it(self.P_sat),
            v_l=it(self.v_l),  v_g=it(self.v_g),
            u_l=it(self.u_l),  u_g=it(self.u_g),
            h_l=it(self.h_l),  h_g=it(self.h_g),
            cp_l=it(self.cp_l), cp_g=it(self.cp_g),
            cv_l=it(self.cv_l), cv_g=it(self.cv_g),
            k_l=it(self.k_l),  k_g=it(self.k_g),
        )


def load_saturation_table(path: str | Path) -> SatProps:
    """Load a NIST saturation table (.cgi, tab-separated). Robust to undefined cells."""
    df = pd.read_csv(path, sep="\t", engine="python")

    T   = _to_num(df["Temperature (K)"])
    P   = _to_num(df["Pressure (MPa)"], factor=1e6)
    v_l = _to_num(df["Volume (l, m3/kg)"])
    v_g = _to_num(df["Volume (v, m3/kg)"])
    u_l = _to_num(df["Internal Energy (l, kJ/kg)"], factor=1e3)
    u_g = _to_num(df["Internal Energy (v, kJ/kg)"], factor=1e3)
    h_l = _to_num(df["Enthalpy (l, kJ/kg)"], factor=1e3)
    h_g = _to_num(df["Enthalpy (v, kJ/kg)"], factor=1e3)
    cp_l = _to_num(df["Cp (l, J/g*K)"], factor=1e3)
    cp_g = _to_num(df["Cp (v, J/g*K)"], factor=1e3)
    cv_l = _to_num(df["Cv (l, J/g*K)"], factor=1e3)
    cv_g = _to_num(df["Cv (v, J/g*K)"], factor=1e3)
    k_l  = _to_num(df["Therm. Cond. (l, W/m*K)"])
    k_g  = _to_num(df["Therm. Cond. (v, W/m*K)"])

    idx = np.argsort(T)
    T = T[idx]
    arrays = (P, v_l, v_g, u_l, u_g, h_l, h_g, cp_l, cp_g, cv_l, cv_g, k_l, k_g)
    P, v_l, v_g, u_l, u_g, h_l, h_g, cp_l, cp_g, cv_l, cv_g, k_l, k_g = (
        a[idx] for a in arrays
    )
    fill = lambda a: _interp_fill_by_T(T, a)
    return SatProps(
        T_grid=T,
        P_sat=fill(P),
        v_l=fill(v_l),  v_g=fill(v_g),
        u_l=fill(u_l),  u_g=fill(u_g),
        h_l=fill(h_l),  h_g=fill(h_g),
        cp_l=fill(cp_l), cp_g=fill(cp_g),
        cv_l=fill(cv_l), cv_g=fill(cv_g),
        k_l=fill(k_l),  k_g=fill(k_g),
    )


# ---------------------------------------------------------------------------
# Isothermal branches
# ---------------------------------------------------------------------------

@dataclass
class IsoBranch:
    """One isotherm restricted to a single phase, with monotone P(rho)."""
    T:     float
    phase: str
    rho:   np.ndarray
    P:     np.ndarray
    u:     np.ndarray
    h:     np.ndarray
    cp:    np.ndarray
    cv:    np.ndarray
    k:     np.ndarray


@dataclass
class IsoAtT:
    """Two phase branches (liquid / gas) at a single temperature."""
    T:   float
    liq: Optional[IsoBranch]
    gas: Optional[IsoBranch]


def load_isothermal_table(path: str | Path) -> pd.DataFrame:
    """Load a single NIST isothermal table; return unified DataFrame in SI."""
    df = pd.read_csv(path, sep="\t", engine="python")

    def pick(prefix: str) -> str:
        for c in df.columns:
            if c.lower().startswith(prefix.lower()):
                return c
        raise KeyError(f"Column starting with '{prefix}' not found in {path}")

    cP = pick("Pressure")
    facP = 1.0
    if "MPa" in cP:   facP = 1e6
    elif "kPa" in cP: facP = 1e3
    elif "atm" in cP: facP = 101_325.0
    P = _to_num(df[cP], factor=facP)

    rho = _to_num(df[pick("Density")])
    cU = pick("Internal Energy")
    cH = pick("Enthalpy")
    u  = _to_num(df[cU], factor=1e3 if "kJ/kg" in cU else 1.0)
    h  = _to_num(df[cH], factor=1e3 if "kJ/kg" in cH else 1.0)
    cCP, cCV = pick("Cp"), pick("Cv")
    cp = _to_num(df[cCP], factor=1e3 if "J/g*K" in cCP else 1.0)
    cv = _to_num(df[cCV], factor=1e3 if "J/g*K" in cCV else 1.0)
    k  = _to_num(df[pick("Therm. Cond.")])

    phase = None
    for c in df.columns:
        if c.strip().lower() == "phase":
            phase = df[c].astype(str).str.lower()
            break

    return pd.DataFrame({
        "P_Pa": P, "rho_kg_m3": rho,
        "u_J_kg": u, "h_J_kg": h,
        "cp_J_kgK": cp, "cv_J_kgK": cv, "k_W_mK": k,
        "Phase": phase if phase is not None else "",
    })


def preprocess_isothermal(df_iso: pd.DataFrame, T: float, sat: SatProps) -> IsoAtT:
    """Split an isotherm into monotone liquid/gas branches, mapping props onto P."""
    df = df_iso.dropna(subset=["P_Pa", "rho_kg_m3"]).copy().sort_values("P_Pa")

    mask_two = (
        df["Phase"].str.contains("two", na=False)
        | df["Phase"].str.contains("sat", na=False)
        | df["Phase"].str.contains(r"l\+v", regex=True, na=False)
    )
    df = df[~mask_two]

    s = sat.eval(T)
    v = 1.0 / df["rho_kg_m3"].to_numpy()
    is_liq = v < s["v_l"]
    is_gas = v > s["v_g"]

    branches: Dict[str, Optional[IsoBranch]] = {}
    for key, mask in [("liq", is_liq), ("gas", is_gas)]:
        if mask.sum() < 3:
            branches[key] = None
            continue
        sub = df.loc[mask]
        rho_arr = sub["rho_kg_m3"].to_numpy()
        P_arr   = sub["P_Pa"].to_numpy()

        rho_m, P_m = _unique_monotone(rho_arr, P_arr, increasing=True)

        def map_to_P_grid(P_src: np.ndarray, Y_src: np.ndarray) -> np.ndarray:
            P_u, Y_u = _unique_monotone(P_src, Y_src, increasing=True)
            return _interp1(P_m, P_u, Y_u)

        branches[key] = IsoBranch(
            T=T, phase=key, rho=rho_m, P=P_m,
            u=map_to_P_grid(P_arr, sub["u_J_kg"].to_numpy()),
            h=map_to_P_grid(P_arr, sub["h_J_kg"].to_numpy()),
            cp=map_to_P_grid(P_arr, sub["cp_J_kgK"].to_numpy()),
            cv=map_to_P_grid(P_arr, sub["cv_J_kgK"].to_numpy()),
            k=map_to_P_grid(P_arr, sub["k_W_mK"].to_numpy()),
        )

    return IsoAtT(T=T, liq=branches["liq"], gas=branches["gas"])


# ---------------------------------------------------------------------------
# PropertyDB
# ---------------------------------------------------------------------------

class PropertyDB:
    """Bundled saturation + isothermal database with state evaluation."""

    def __init__(self, sat: SatProps, iso_map: Dict[float, IsoAtT]):
        if not iso_map:
            raise ValueError("iso_map is empty; provide at least one isothermal table")
        self.sat = sat
        self.T_list  = sorted(iso_map.keys())
        self.iso_at_T = iso_map

    # ------ construction helpers ------

    @staticmethod
    def _parse_T_from_path(path: str | Path) -> Optional[float]:
        m = re.search(r"(\d+(?:\.\d+)?)\s*K", os.path.basename(str(path)))
        return float(m.group(1)) if m else None

    @classmethod
    def from_files(cls, sat_path: str | Path,
                   isotherm_paths: Iterable[str | Path]) -> "PropertyDB":
        """Build a PropertyDB from explicit file paths."""
        sat = load_saturation_table(sat_path)
        iso_map: Dict[float, IsoAtT] = {}
        for p in isotherm_paths:
            df_iso = load_isothermal_table(p)
            T_val = cls._parse_T_from_path(p)
            if T_val is None:
                raise ValueError(f"Cannot infer temperature from filename: {p}")
            iso_map[float(T_val)] = preprocess_isothermal(df_iso, float(T_val), sat)
        return cls(sat, iso_map)

    @classmethod
    def from_directory(cls, sat_path: str | Path, isotherm_dir: str | Path
                       ) -> "PropertyDB":
        """Build a PropertyDB from a directory of isothermal .cgi files."""
        files = sorted(Path(isotherm_dir).glob("*.cgi"))
        if not files:
            raise FileNotFoundError(f"No .cgi files found under {isotherm_dir}")
        return cls.from_files(sat_path, files)

    # ------ state evaluation ------

    def _bracket_T(self, T: float) -> Tuple[float, float]:
        if T <= self.T_list[0]:  return self.T_list[0], self.T_list[0]
        if T >= self.T_list[-1]: return self.T_list[-1], self.T_list[-1]
        for i in range(len(self.T_list) - 1):
            if self.T_list[i] <= T <= self.T_list[i + 1]:
                return self.T_list[i], self.T_list[i + 1]
        return self.T_list[-1], self.T_list[-1]

    def state_from_T_rho(self, T: float, rho: float) -> Dict[str, float]:
        """Evaluate the equation of state at (T, rho).

        Returns a dictionary with keys
        ``T, rho, P, u, h, cp, cv, k, gamma, phase`` (and ``x`` for two-phase).
        """
        T = float(T); rho = float(rho)
        satT = self.sat.eval(T)
        v    = 1.0 / rho
        v_l  = satT["v_l"]
        v_g  = satT["v_g"]

        # Two-phase region: clamp P to P_sat, lever-rule for extensive props
        if v_l < v < v_g:
            x = (v - v_l) / (v_g - v_l)
            blend = lambda a, b: (1.0 - x) * a + x * b
            cp = blend(satT["cp_l"], satT["cp_g"])
            cv = blend(satT["cv_l"], satT["cv_g"])
            return dict(
                T=T, rho=rho, P=satT["P_sat"],
                u=blend(satT["u_l"], satT["u_g"]),
                h=blend(satT["h_l"], satT["h_g"]),
                cp=cp, cv=cv, k=blend(satT["k_l"], satT["k_g"]),
                gamma=cp / max(cv, 1e-9),
                phase="two-phase", x=x,
            )

        # Single-phase: bilinear interpolation on the isothermal tables
        T0, T1 = self._bracket_T(T)
        w = 0.0 if T0 == T1 else (T - T0) / (T1 - T0)

        def eval_single(atT: Optional[IsoAtT]) -> Optional[Dict[str, float]]:
            if atT is None:
                return None
            sat_loc = self.sat.eval(atT.T)
            branch = (atT.liq or atT.gas) if v <= sat_loc["v_l"] else (atT.gas or atT.liq)
            if branch is None:
                return None
            rhoq = _clamp(rho, float(branch.rho[0]), float(branch.rho[-1]))
            Pq = float(_interp1(np.array([rhoq]), branch.rho, branch.P)[0])
            atP = lambda arr: float(_interp1(np.array([Pq]), branch.P, arr)[0])
            u, h, cp, cv, k = (atP(branch.u), atP(branch.h),
                               atP(branch.cp), atP(branch.cv), atP(branch.k))
            return dict(P=Pq, u=u, h=h, cp=cp, cv=cv, k=k,
                        gamma=cp / max(cv, 1e-9), phase=branch.phase)

        s0 = eval_single(self.iso_at_T.get(T0))
        s1 = eval_single(self.iso_at_T.get(T1))

        if s0 is None and s1 is None:
            # Fallback: use saturation gas at nearest T
            s_near = self.sat.eval(T0)
            cp = s_near["cp_g"]
            cv = s_near["cv_g"]
            return dict(T=T, rho=rho, P=s_near["P_sat"],
                        u=s_near["u_g"], h=s_near["h_g"],
                        cp=cp, cv=cv, k=s_near["k_g"],
                        gamma=cp / max(cv, 1e-9),
                        phase="fallback", x=float("nan"))

        if s0 is None: s0 = s1
        if s1 is None: s1 = s0

        out = {key: (1.0 - w) * s0[key] + w * s1[key]
               for key in ["P", "u", "h", "cp", "cv", "k", "gamma"]}
        out.update(T=T, rho=rho, x=float("nan"),
                   phase="liq" if v <= v_l else "gas")
        return out


__all__ = [
    "SatProps", "IsoBranch", "IsoAtT",
    "PropertyDB",
    "load_saturation_table", "load_isothermal_table", "preprocess_isothermal",
]
