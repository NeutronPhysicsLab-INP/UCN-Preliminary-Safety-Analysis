"""
Loss-of-Vacuum Accident (LOVA) with helium coolant gap (Section 3.6).

Air ingress into the vacuum jacket triggers cryopumping deposition on the
cold surfaces; the heat propagates through the helium-jacket outer wall,
across the trapped helium gap (which pressurizes isochorically), through
the outer deuterium wall, and into the LD2 inventory. The rupture disk
actuates when the deuterium pressure reaches its setpoint.

Two coupled lumped-parameter wall-temperature ODEs (outer wall + inner
wall) are integrated together with the deuterium energy balance and the
ballast-tank dynamics.
"""

from __future__ import annotations
import math
from typing import Dict, List, Tuple

import pandas as pd

from ..config import SimConfig
from ..eos.solver import find_rho_for_target_P
from ..relief.rupture_disk import step_burst_disk_absolute


def get_cp_aluminum(T: float) -> float:
    """Specific heat of Aluminum 5056/6061 (approximate, 4 - 300 K).

    Curve fit used in the manuscript LOVA model. Returns J/(kg K).
    """
    T = max(4.0, T)
    if T < 40.0:  return 0.9 * T
    if T < 100.0: return 20.0 + 6.0 * (T - 40.0)
    if T < 300.0: return 480.0 + 2.1 * (T - 100.0)
    return 900.0


def _chamber_areas(geom) -> Tuple[float, float, float, float]:
    """Return (A_cyl_outer, A_sph_outer, gap_cyl, gap_sph) for the LOVA load."""
    R_outer_face = geom.R_he_jacket + geom.t_he_jacket
    L_cyl = geom.L_cyl
    angle_rad = math.radians(geom.cap_angle_deg)

    # Cylindrical body outer area
    A_cyl_outer = 2.0 * math.pi * R_outer_face * L_cyl

    # Spherical-cap outer area
    factor = (1.0 - math.cos(angle_rad))
    A_sph_outer = 2.0 * math.pi * (R_outer_face ** 2) * factor

    # Vacuum gaps (cylindrical and spherical)
    gap_cyl = geom.R_vacuum_inner - R_outer_face
    gap_sph = geom.R_vacuum_inner - R_outer_face  # same convention as notebook

    return A_cyl_outer, A_sph_outer, gap_cyl, gap_sph


def calc_lova_heat_load(t: float, T_wall: float, *,
                        A_cyl_outer: float, A_sph_outer: float,
                        gap_cyl: float, gap_sph: float,
                        T_outer: float = 298.0,
                        flux_dep: float = 38_000.0,
                        eps_rad:  float = 0.055,
                        P_half_knudsen: float = 100.0) -> Tuple[float, float]:
    """External heat load Q_ext(t, T_wall) for the LOVA scenario.

    Sums residual-gas conduction (Knudsen-flow interpolation), natural
    convection, cryopumping deposition, and Stefan-Boltzmann radiation.
    """
    # Vacuum ramp 1.3e-3 -> 101.3 kPa over 100 s (slope hard-coded for now;
    # could be exposed via cfg.lova.pressurization_t)
    P_vac = min(1.3e-3 + 1000.0 * t, 101_325.0)
    k_atm = 0.025
    k_eff = k_atm * P_vac / (P_vac + P_half_knudsen)

    # 1. Residual gas conduction
    Q_cond_cyl = k_eff * A_cyl_outer * (T_outer - T_wall) / gap_cyl
    Q_cond_sph = k_eff * A_sph_outer * (T_outer - T_wall) / gap_sph

    # 2. Natural convection
    Q_conv_cyl = Q_conv_sph = 0.0
    if P_vac > 500.0:
        deltaT = max(0.0, T_outer - T_wall)
        h_conv = 5.0 + 0.05 * (P_vac - 500.0) ** 0.5
        Q_conv_cyl = h_conv * A_cyl_outer * deltaT
        Q_conv_sph = h_conv * A_sph_outer * deltaT

    # 3. Cryopumping deposition
    Q_dep = 0.0
    if P_vac > 10.0:
        p_scale = min(1.0, (P_vac - 10.0) / 10_000.0)
        if T_wall < 65.0:
            t_scale = 1.0
        elif T_wall > 80.0:
            t_scale = 0.0
        else:
            t_scale = (80.0 - T_wall) / 15.0
        Q_dep = flux_dep * (A_cyl_outer + A_sph_outer) * p_scale * t_scale

    # 4. Stefan-Boltzmann radiation (polished cryogenic Al)
    sigma_sb = 5.67e-8
    Q_rad = sigma_sb * (A_cyl_outer + A_sph_outer) * (T_outer ** 4 - T_wall ** 4) * eps_rad

    return Q_cond_cyl + Q_conv_cyl + Q_cond_sph + Q_conv_sph + Q_dep + Q_rad, P_vac


def run_lova_with_helium_gap(
    db,
    cfg: SimConfig,
    *,
    m1_0: float = None,
    T1_0: float = None,
    P_burst_abs: float = None,
    Q_scale: float = 1.0,
    t_end: float = None,
    mass_inner: float,
    mass_outer: float,
) -> pd.DataFrame:
    """LOVA simulation with the helium-jacket coolant gap.

    Solves coupled ODEs for the outer wall, inner wall, helium-gap pressure,
    and the deuterium-vessel + ballast-tank dynamics.

    Parameters
    ----------
    mass_inner, mass_outer : aluminium-shell masses (kg) for the wall
        thermal-inertia ODEs. Compute these once from the geometry and
        density of Al 5056 (see ``scripts/run_lova_with_helium.py``).

    Returns
    -------
    pandas DataFrame with columns including ``t, P1, P_he, T_in, T_out,
    m1, m2, T1, T2, Q_ext`` (sampled every ~5 s, plus dense sampling
    across the disk-burst event).
    """
    geom  = cfg.geometry
    deut  = cfg.deuterium
    ruptd = cfg.rupture_disk
    he    = cfg.helium

    if m1_0 is None:        m1_0 = deut.mass
    if T1_0 is None:        T1_0 = deut.T_init
    if P_burst_abs is None: P_burst_abs = ruptd.p_burst
    if t_end is None:       t_end = cfg.uq.sim_time_lova

    V1 = geom.V_inventory
    V2 = ruptd.V_ballast

    A_cyl_outer, A_sph_outer, gap_cyl, gap_sph = _chamber_areas(geom)

    # Init fluids
    rho1 = m1_0 / V1
    s1 = db.state_from_T_rho(T1_0, rho1)
    m1, U1, T1, P1 = m1_0, s1["u"] * m1_0, T1_0, s1["P"]

    # Walls (start at deuterium temperature)
    T_w_in, T_w_out = T1_0, T1_0

    # Trapped helium gap: fixed mass at initial conditions
    P_he_init = he.p_init
    R_spec_he = he.R_specific
    rho_he_gap = P_he_init / (R_spec_he * T1_0)

    # Ballast tank
    rho2 = find_rho_for_target_P(db, ruptd.T_ballast_init, ruptd.p_ballast_init)
    m2 = rho2 * V2
    s2 = db.state_from_T_rho(ruptd.T_ballast_init, rho2)
    U2, T2 = s2["u"] * m2, ruptd.T_ballast_init

    t, dt = 0.0, 0.001
    valve_open, time_since = False, 0.0

    # Conductances (constants in the lumped model)
    h_boil = 5_000.0
    A_inner_total = 0.50  # m^2 — approximate inner area
    R_cond_wall = geom.t_outer_wall / (237.0 * A_inner_total)  # k_Al ≈ 237 W/m/K at low T
    R_conv_film = 1.0 / (h_boil * A_inner_total)

    hist: List[Dict] = [
        dict(t=0.0, P1=P1, P_he=P_he_init,
             T_in=T_w_in, T_out=T_w_out, Q_ext=0.0)
    ]

    while t < t_end:
        # External heat onto outer wall
        Q_ext_base, _ = calc_lova_heat_load(
            t, T_w_out,
            A_cyl_outer=A_cyl_outer, A_sph_outer=A_sph_outer,
            gap_cyl=gap_cyl, gap_sph=gap_sph,
            flux_dep=cfg.lova.flux_deposition,
            eps_rad=cfg.lova.epsilon_radiative,
            P_half_knudsen=cfg.lova.P_half_knudsen,
        )
        Q_ext = Q_ext_base * Q_scale

        # Helium-gap heat transfer
        T_he_avg = 0.5 * (T_w_out + T_w_in)
        k_he_local = 0.02 + 0.0004 * T_he_avg
        R_he = geom.he_gap / (k_he_local * A_inner_total)
        Q_he = (T_w_out - T_w_in) / R_he

        # Helium-gap pressure (isochoric ideal-gas)
        P_he = rho_he_gap * R_spec_he * T_he_avg

        # Heat into the deuterium fluid
        Q_fluid = (T_w_in - T1) / (R_cond_wall + R_conv_film)

        # Update walls
        cp_out = get_cp_aluminum(T_w_out)
        cp_in  = get_cp_aluminum(T_w_in)
        T_w_out += (Q_ext - Q_he) * dt / (mass_outer * cp_out)
        T_w_in  += (Q_he - Q_fluid) * dt / (mass_inner * cp_in)

        # Burst-disk step
        step = step_burst_disk_absolute(
            db, m1, T1, U1, m2, T2, U2, t, dt, V1, V2,
            mdot_open=15.0,
            P_burst_limit=P_burst_abs,
            Q1_in=Q_fluid, Q2_loss=0.0,
            diameter_m=ruptd.diameter,
            Cd=ruptd.Cd, t_open=ruptd.open_time,
            valve_open=valve_open, time_since_burst=time_since,
        )

        m1, T1, U1 = step["m1"], step["T1"], step["U1"]
        m2, T2, U2 = step["m2"], step["T2"], step["U2"]
        valve_open  = step["valve_open"]
        time_since  = step["time_since_open"]
        t           = step["t"]

        # Adaptive time-stepping
        if valve_open:
            dt = 0.0005 if time_since < 2.0 else 0.02
        else:
            if abs((Q_ext - Q_he) * dt / (mass_outer * cp_out)) > 2.0:
                dt *= 0.5
            elif dt < 0.01:
                dt *= 1.1

        # Sample roughly every 5 s, densely across burst
        if t % 5 < dt or (valve_open and time_since < 5.0):
            step.update(T_out=T_w_out, T_in=T_w_in, Q_ext=Q_ext, P_he=P_he)
            hist.append(step)

    return pd.DataFrame(hist)


__all__ = [
    "calc_lova_heat_load", "run_lova_with_helium_gap", "get_cp_aluminum",
]
