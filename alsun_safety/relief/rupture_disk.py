"""
Rupture-disk dynamics for the deuterium-vessel / ballast-tank relief path.

The disk is modeled as a discrete-event trigger:

    1. While ``P_vessel < P_set`` it remains closed (no flow).
    2. The first time ``P_vessel >= P_set`` it bursts; the flow area then
       ramps from 0 to A_full over a characteristic opening time
       ``t_open`` (default 0.1 s).
    3. After full opening, mass flow is governed by the Bernoulli
       (subcooled liquid) or API 520 compressible-flow equations.

The single time-step routine :func:`step_burst_disk_absolute` advances
one explicit-Euler step of the lumped-parameter conservation laws across
the rupture-disk interface.
"""

from __future__ import annotations
from typing import Callable, Dict, Union

import numpy as np

from .flow import _as_func, mdot_bernoulli, mdot_orifice
from ..eos.solver import invert_T_from_u


def step_burst_disk_absolute(
    db,
    m1: float, T1: float, U1: float,
    m2: float, T2: float, U2: float,
    t: float, dt: float,
    V1: float, V2: float,
    *,
    mdot_open: float,
    P_burst_limit: float,
    Q1_in: Union[float, Callable[[float], float]],
    Q2_loss: Union[float, Callable[[float], float]],
    diameter_m: float,
    Cd: float = 0.8,
    t_open: float = 0.1,
    efflux_model: str = "vapor_only",
    valve_open: bool = False,
    time_since_burst: float = 0.0,
) -> Dict[str, float]:
    """Advance the deuterium vessel + ballast tank by one time step.

    Parameters
    ----------
    db                : :class:`PropertyDB` for the working fluid.
    m1, T1, U1        : vessel-1 (deuterium vessel) state at time t.
    m2, T2, U2        : vessel-2 (ballast tank) state at time t.
    V1, V2            : vessel volumes (m^3).
    mdot_open         : pipe-choke ceiling on mass-flow rate (kg/s).
    P_burst_limit     : disk activation pressure (Pa, absolute).
    Q1_in, Q2_loss    : heat input (W) to vessel 1, heat loss (W) from
                        vessel 2; either constants or functions of t.
    diameter_m        : disk diameter (m). 25 mm = 0.025.
    Cd                : discharge coefficient.
    t_open            : characteristic opening time (s).
    efflux_model      : ``"vapor_only"`` clamps two-phase efflux to vapor
                        properties (used in the manuscript); anything else
                        uses the local two-phase mixture.
    valve_open        : has the disk already burst on a prior step?
    time_since_burst  : seconds since burst (controls area ramp).

    Returns
    -------
    Dictionary describing the new state at time t+dt.
    """
    Q1f = _as_func(Q1_in)
    Q2f = _as_func(Q2_loss)

    rho1 = m1 / V1
    s1 = db.state_from_T_rho(T1, rho1)
    P1 = s1["P"]

    rho2 = m2 / V2
    s2 = db.state_from_T_rho(T2, rho2)
    P2 = s2["P"]

    dP = P1 - P2

    # Burst trigger
    if not valve_open:
        if P1 >= P_burst_limit:
            valve_open = True
            time_since_burst = 0.0
    else:
        time_since_burst += dt

    mdot = 0.0
    h_out = s1["h"]

    if valve_open and dP > 0:
        # Linear area ramp from 0 to full over t_open seconds
        s_val = min(1.0, time_since_burst / t_open)
        A_eff = (np.pi * (diameter_m / 2.0) ** 2) * s_val

        if s1["phase"] in ("liq", "liq_sat_clamped"):
            mdot_phys = mdot_bernoulli(rho1, P1, P2, A_eff, Cd)
            h_out = s1["h"]
        else:
            if efflux_model == "vapor_only" and s1["phase"] == "two-phase":
                sat = db.sat.eval(T1)
                cp_f, cv_f, h_out = sat["cp_g"], sat["cv_g"], sat["h_g"]
            else:
                cp_f, cv_f, h_out = s1["cp"], s1["cv"], s1["h"]
            mdot_phys = mdot_orifice(P1, T1, P2, A_eff, Cd, cp_f, cv_f)

        mdot = min(mdot_phys, mdot_open)

    # Mass / energy update
    dm = mdot * dt
    if dm > m1:
        dm, mdot = m1, m1 / dt

    m1n = m1 - dm
    m2n = m2 + dm

    # Cut external heating once the vessel is essentially empty so the EOS
    # solver does not chase a singular high-pressure state.
    Q_eff = 0.0 if m1n < 1e-4 else Q1f(t)

    U1n = U1 + (Q_eff - h_out * mdot) * dt
    U2n = U2 + (h_out * mdot - Q2f(t)) * dt

    rho1n = m1n / V1
    if m1n > 1e-4:
        T1n = invert_T_from_u(db, rho1n, U1n / m1n, T_guess=T1)
    else:
        T1n = T2  # equilibrium with ballast

    rho2n = m2n / V2
    T2n = invert_T_from_u(db, rho2n, U2n / m2n, T_guess=T2)

    s1n = db.state_from_T_rho(T1n, rho1n)
    s2n = db.state_from_T_rho(T2n, rho2n)

    return dict(
        t=t + dt,
        m1=m1n, T1=T1n, U1=U1n, P1=s1n["P"], phase1=s1["phase"],
        m2=m2n, T2=T2n, U2=U2n, P2=s2n["P"],
        mdot=mdot,
        valve_open=valve_open,
        time_since_open=time_since_burst,
    )


__all__ = ["step_burst_disk_absolute"]
