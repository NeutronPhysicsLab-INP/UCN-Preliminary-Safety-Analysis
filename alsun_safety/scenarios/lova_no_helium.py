"""
Loss-of-Vacuum Accident (LOVA) without the helium coolant gap — single-shell
variant used as a counterfactual comparison.

The deuterium chamber is treated as a single-walled vessel directly in
contact with the vacuum environment (i.e., the helium-jacket buffer is
absent). The full air-ingress / cryopumping heat load lands on the single
deuterium-bounding wall, and the deuterium fluid sees the heat after
thermal-inertia delay through that wall and a nucleate-boiling resistance.
This represents alternative designs with a circulating deuterium loop
connected directly to an external heat exchanger, where no helium buffer
is present.

Only one simulator function is exposed: :func:`run_lova_without_helium`.
The natural-convection / conduction / cryopumping / radiation breakdown
re-uses :func:`alsun_safety.scenarios.lova_helium.calc_lova_heat_load` so
that the ``with-He`` and ``without-He`` scenarios share the same external
loads.
"""

from __future__ import annotations
import math
from typing import Dict, List

import pandas as pd

from ..config import SimConfig
from ..eos.solver import find_rho_for_target_P
from ..relief.rupture_disk import step_burst_disk_absolute
from .lova_helium import calc_lova_heat_load, get_cp_aluminum, _chamber_areas


def run_lova_without_helium_gap(
    db,
    cfg: SimConfig,
    *,
    m1_0: float = None,
    T1_0: float = None,
    P_burst_abs: float = None,
    Q_scale: float = 1.0,
    t_end: float = None,
    mass_wall: float = None,
) -> pd.DataFrame:
    """LOVA simulation for a single-walled vessel (no helium coolant gap).

    Parameters
    ----------
    mass_wall : float, optional
        Aluminium-shell mass (kg) for the single-wall thermal-inertia ODE.
        If ``None``, computed from cfg.geometry assuming the
        ``inner deuterium chamber wall`` thickness is the load-bearing
        single shell.

    Returns
    -------
    pandas DataFrame with columns ``t, P1, T1, T_wall, m1, m2, T2,
    Q_in, Q_ext`` (sampled every ~5 s plus densely across the disk-burst
    event).
    """
    geom = cfg.geometry
    deut = cfg.deuterium
    ruptd = cfg.rupture_disk

    if m1_0 is None:        m1_0 = deut.mass
    if T1_0 is None:        T1_0 = deut.T_init
    if P_burst_abs is None: P_burst_abs = ruptd.p_burst
    if t_end is None:       t_end = cfg.uq.sim_time_lova

    V1 = geom.V_inventory
    V2 = ruptd.V_ballast

    # Compute single-shell wall mass if not given. Use the inner deuterium
    # chamber wall thickness as the relevant shell.
    if mass_wall is None:
        rho_al = cfg.aluminum.rho
        L_cyl  = geom.L_cyl
        t_wall = geom.t_inner_wall
        R      = geom.R_inner_wall
        # Cylindrical body + (thin) hemispherical-cap approximation
        vol_cyl = math.pi * ((R + t_wall) ** 2 - R ** 2) * L_cyl
        vol_cap = 4.0 * math.pi * R ** 2 * t_wall / 2.0  # ~hemisphere shell
        mass_wall = (vol_cyl + vol_cap) * rho_al

    A_cyl_outer, A_sph_outer, gap_cyl, gap_sph = _chamber_areas(geom)

    # Initial state
    rho1 = m1_0 / V1
    s1 = db.state_from_T_rho(T1_0, rho1)
    m1, U1, T1, P1 = m1_0, s1["u"] * m1_0, T1_0, s1["P"]
    T_wall = T1_0

    # Ballast tank
    rho2 = find_rho_for_target_P(db, ruptd.T_ballast_init, ruptd.p_ballast_init)
    m2 = rho2 * V2
    s2 = db.state_from_T_rho(ruptd.T_ballast_init, rho2)
    U2, T2 = s2["u"] * m2, ruptd.T_ballast_init

    t, dt = 0.0, 0.01
    valve_open, time_since = False, 0.0

    # Conductance through the single wall + boiling film
    h_boil = 5_000.0
    A_inner = 0.50  # m^2, approximate
    R_cond = geom.t_inner_wall / (237.0 * A_inner)
    R_conv = 1.0 / (h_boil * A_inner)
    R_total = R_cond + R_conv

    hist: List[Dict] = [dict(t=0.0, P1=P1, T1=T1, T_wall=T_wall, Q_in=0.0, Q_ext=0.0)]

    while t < t_end:
        Q_ext_base, _ = calc_lova_heat_load(
            t, T_wall,
            A_cyl_outer=A_cyl_outer, A_sph_outer=A_sph_outer,
            gap_cyl=gap_cyl, gap_sph=gap_sph,
            flux_dep=cfg.lova.flux_deposition,
            eps_rad=cfg.lova.epsilon_radiative,
            P_half_knudsen=cfg.lova.P_half_knudsen,
        )
        Q_ext = Q_ext_base * Q_scale

        # Heat into the deuterium fluid
        Q_fluid = (T_wall - T1) / R_total

        # Update wall temperature (single-shell thermal inertia)
        cp_al = get_cp_aluminum(T_wall)
        T_wall += (Q_ext - Q_fluid) * dt / (mass_wall * cp_al)

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
            dt = 0.001 if time_since < 2.0 else 0.05

        if t % 5 < dt or (valve_open and time_since < 5.0):
            step.update(T_wall=T_wall, Q_in=Q_fluid, Q_ext=Q_ext)
            hist.append(step)

    return pd.DataFrame(hist)


__all__ = ["run_lova_without_helium_gap"]
