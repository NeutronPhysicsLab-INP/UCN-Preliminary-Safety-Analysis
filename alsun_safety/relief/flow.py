"""
Orifice flow models used by the rupture-disk module.

Implements:

- Bernoulli (incompressible) flow for the subcooled-liquid phase of the
  blowdown (Equation in §3.2 of the manuscript).
- API Standard 520 (Part I) compressible-flow formulation, automatically
  switching between choked (sonic) and subsonic regimes via the critical
  pressure ratio :math:`r_{crit} = \\big( 2/(\\gamma+1) \\big)^{\\gamma/(\\gamma-1)}`.
"""

from __future__ import annotations
import numpy as np


def _as_func(x):
    """Wrap a constant scalar as a callable f(t)->float, or return the callable."""
    return x if callable(x) else (lambda t: float(x))


def mdot_bernoulli(rho_up: float, P_up: float, P_down: float,
                   A: float, Cd: float) -> float:
    """Incompressible Bernoulli mass-flow rate (subcooled-liquid blowdown phase).

    .. math::
        \\dot{m}_{liq} = C_d A \\sqrt{2 \\rho_{up} (P_{up} - P_{down})}
    """
    if A <= 0.0 or P_up <= P_down:
        return 0.0
    return Cd * A * float(np.sqrt(2.0 * rho_up * (P_up - P_down)))


def mdot_orifice(P_up: float, T_up: float, P_down: float,
                 A: float, Cd: float, cp: float, cv: float) -> float:
    """Compressible-gas orifice flow (API Standard 520, Part I).

    Switches automatically between choked (sonic) and subsonic flow.

    Parameters
    ----------
    P_up, T_up : upstream pressure (Pa) and temperature (K).
    P_down     : downstream pressure (Pa).
    A          : effective flow area (m^2).
    Cd         : discharge coefficient.
    cp, cv     : specific heats (J / kg / K) -> gamma, R inferred.

    Returns
    -------
    Mass flow rate (kg/s); zero if flow direction is reversed.
    """
    if A <= 0.0 or P_up <= P_down:
        return 0.0
    gamma = cp / max(cv, 1e-9)
    R = max(cp - cv, 1e-9)
    pr = P_down / P_up
    pr_crit = (2.0 / (gamma + 1.0)) ** (gamma / (gamma - 1.0))

    if pr <= pr_crit:
        # Choked (sonic)
        factor = (2.0 / (gamma + 1.0)) ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))
        m_dot = Cd * A * P_up * np.sqrt(gamma / (R * T_up)) * factor
    else:
        # Subsonic
        term = pr ** (2.0 / gamma) - pr ** ((gamma + 1.0) / gamma)
        term = max(term, 0.0)
        m_dot = Cd * A * P_up * np.sqrt(2.0 * gamma / (R * T_up * (gamma - 1.0)) * term)

    return float(max(m_dot, 0.0))


__all__ = ["_as_func", "mdot_bernoulli", "mdot_orifice"]
