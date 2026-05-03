"""
Cooling-failure scenario simulator (Section 3.5 of the manuscript).

Loss of helium refrigeration with a constant heat-load source ``Q``
applied to the deuterium inventory. The deuterium pressure rises until
the rupture disk actuates at ``P_burst_abs``, blowing down into a warm
ballast tank initialized at atmospheric pressure.

Returns a dict-of-lists time history of ``t, P1, P2, m1, m2, T1, T2``
that downstream stress and plotting routines accept directly.
"""

from __future__ import annotations
from typing import Callable, Dict, List, Union

from ..config import SimConfig
from ..eos.solver import find_rho_for_target_P
from ..relief.rupture_disk import step_burst_disk_absolute


def run_cooling_failure(
    db,
    cfg: SimConfig,
    *,
    m1_0: float = None,
    T1_0: float = None,
    Q_total: Union[float, Callable[[float], float]] = None,
    P_burst_abs: float = None,
    t_end: float = None,
) -> Dict[str, List[float]]:
    """Run the cooling-failure transient.

    Parameters
    ----------
    db : :class:`PropertyDB`
        Loaded equation-of-state database.
    cfg : :class:`SimConfig`
        Manuscript-baseline configuration. Most arguments default to
        the corresponding fields of ``cfg``.
    m1_0, T1_0 : initial deuterium mass (kg) and temperature (K).
                 Defaults: cfg.deuterium.mass, cfg.deuterium.T_init.
    Q_total    : heat load (W) on the deuterium inventory (constant or
                 callable f(t)). Default: cfg.cooling_failure.Q_total.
    P_burst_abs: rupture-disk activation pressure (Pa, absolute).
                 Default: cfg.rupture_disk.p_burst.
    t_end      : simulation duration (s). Default: cfg.uq.sim_time.

    Returns
    -------
    Dict[str, List[float]] with keys ``t, P1, P2, m1, m2, T1, T2``.
    Pressures are returned in **bar** for plotting convenience; mass in kg,
    temperature in K, time in s.
    """
    geom = cfg.geometry
    deut = cfg.deuterium
    ruptd = cfg.rupture_disk
    solv = cfg.solver

    if m1_0 is None:        m1_0 = deut.mass
    if T1_0 is None:        T1_0 = deut.T_init
    if Q_total is None:     Q_total = cfg.cooling_failure.Q_total
    if P_burst_abs is None: P_burst_abs = ruptd.p_burst
    if t_end is None:       t_end = cfg.uq.sim_time

    V1 = geom.V_inventory
    V2 = ruptd.V_ballast

    # Initial state for vessel 1 (deuterium)
    rho1 = m1_0 / V1
    s1 = db.state_from_T_rho(T1_0, rho1)
    m1, U1, T1, P1 = m1_0, s1["u"] * m1_0, T1_0, s1["P"]

    # Initial state for vessel 2 (warm ballast tank, near 1 atm at 300 K)
    rho2 = find_rho_for_target_P(db, ruptd.T_ballast_init, ruptd.p_ballast_init)
    m2 = rho2 * V2
    s2 = db.state_from_T_rho(ruptd.T_ballast_init, rho2)
    U2, T2 = s2["u"] * m2, ruptd.T_ballast_init

    t, dt = 0.0, solv.dt_initial
    valve_open, time_since = False, 0.0

    history: Dict[str, List[float]] = {
        "t": [t], "P1": [P1 / 1e5], "P2": [s2["P"] / 1e5],
        "m1": [m1], "m2": [m2], "T1": [T1], "T2": [T2],
    }

    while t < t_end:
        # Adaptive step before burst
        if not valve_open:
            P_dist = P_burst_abs - P1
            if   0 < P_dist < 0.1e5: dt = 0.0005
            elif 0 < P_dist < 0.5e5: dt = 0.002
            else:                    dt = solv.dt_initial

        step = step_burst_disk_absolute(
            db, m1, T1, U1, m2, T2, U2, t, dt, V1, V2,
            mdot_open=ruptd.mdot_max_cap,
            P_burst_limit=P_burst_abs,
            Q1_in=Q_total, Q2_loss=0.0,
            diameter_m=ruptd.diameter,
            Cd=ruptd.Cd, t_open=ruptd.open_time,
            valve_open=valve_open,
            time_since_burst=time_since,
        )

        m1, T1, U1 = step["m1"], step["T1"], step["U1"]
        m2, T2, U2 = step["m2"], step["T2"], step["U2"]
        valve_open, time_since, t, P1 = (
            step["valve_open"], step["time_since_open"], step["t"], step["P1"],
        )

        history["t"].append(t)
        history["P1"].append(step["P1"] / 1e5)
        history["P2"].append(step["P2"] / 1e5)
        history["m1"].append(step["m1"])
        history["m2"].append(step["m2"])
        history["T1"].append(step["T1"])
        history["T2"].append(step["T2"])

        # Adaptive step after burst
        if valve_open:
            if   time_since < 2.0:  dt = 0.001
            elif time_since < 20.0: dt = solv.dt_quasi
            else:                   dt = 1.0

    return history


__all__ = ["run_cooling_failure"]
