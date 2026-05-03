"""
Inversion utilities for the equation-of-state database.

The forward EOS provides P, u, h, ... at a given (T, rho). The transient
solver needs the inverse: given (rho, u_target) find T (used in the energy
update of the lumped-parameter integrator), and given (T, P_target) find
rho (used to set initial conditions).
"""

from __future__ import annotations
from typing import Optional

import numpy as np


def invert_T_from_u(db, rho: float, u_target: float,
                    T_guess: Optional[float] = None,
                    T_min: Optional[float] = None,
                    T_max: Optional[float] = None,
                    tol: float = 1e-3,
                    max_iter: int = 20) -> float:
    """Solve u_EOS(T, rho) = u_target for T.

    Hybrid secant + bisection. Bounds default to the database's loaded
    isotherm range, so temperatures above the critical temperature are
    permitted as long as isotherms exist.
    """
    if T_min is None: T_min = db.T_list[0]
    if T_max is None: T_max = db.T_list[-1]

    if T_guess is None: T_guess = 0.5 * (T_min + T_max)
    T_curr = float(np.clip(T_guess, T_min, T_max))

    def err_at(T_val: float) -> float:
        return db.state_from_T_rho(float(T_val), float(rho))["u"] - u_target

    err_curr = err_at(T_curr)

    # 1) Secant phase
    T_prev   = T_curr * 0.99 if T_curr > T_min else T_curr + 0.5
    err_prev = err_at(T_prev)
    for _ in range(max_iter):
        if abs(err_curr) < tol:
            return T_curr
        denom = err_curr - err_prev
        if abs(denom) < 1e-9:
            break
        delta  = -err_curr * (T_curr - T_prev) / denom
        delta  = max(-10.0, min(10.0, delta))
        T_next = T_curr + delta
        if T_next < T_min: T_next = T_min + 0.1
        if T_next > T_max: T_next = T_max - 0.1
        T_prev, err_prev = T_curr, err_curr
        T_curr = T_next
        err_curr = err_at(T_curr)

    # 2) Bisection fallback
    lo, hi = T_min, T_max
    for _ in range(30):
        mid  = 0.5 * (lo + hi)
        fmid = err_at(mid)
        if abs(fmid) < tol:
            return mid
        # u increases with T at fixed rho
        if fmid > 0: hi = mid
        else:        lo = mid
    return 0.5 * (lo + hi)


def find_rho_for_target_P(db, T_target: float, P_target_Pa: float,
                          *, tol_Pa: float = 100.0,
                          verbose: bool = False) -> float:
    """Solve P_EOS(T, rho) = P_target for rho at fixed T.

    Bisection with a robust upward bracket search. Used to seed initial
    conditions (e.g., rho such that the LD2 vessel sits at exactly 1.5 atm
    at the start of a transient).
    """
    sat = db.sat.eval(T_target)
    rho_sat_l = 1.0 / sat["v_l"]
    P_sat = sat["P_sat"]

    if P_target_Pa <= P_sat:
        if verbose:
            print(f"Target P ({P_target_Pa/1e5:.3f} bar) is below saturation "
                  f"({P_sat/1e5:.3f} bar); returning saturated-liquid density.")
        return rho_sat_l

    # Bracket the target pressure from above
    rho_min = rho_sat_l
    rho_max = rho_sat_l * 1.01
    for i in range(20):
        if db.state_from_T_rho(T_target, rho_max)["P"] > P_target_Pa:
            break
        rho_min = rho_max
        rho_max = rho_max * 1.01
    else:
        if verbose:
            print("Warning: could not bracket target pressure within 20 expansions.")
        return rho_max

    # Bisection
    rho_mid = 0.5 * (rho_min + rho_max)
    P_mid = float("nan")
    for _ in range(50):
        rho_mid = 0.5 * (rho_min + rho_max)
        P_mid = db.state_from_T_rho(T_target, rho_mid)["P"]
        if abs(P_mid - P_target_Pa) < tol_Pa:
            if verbose:
                print(f"Solver converged: P={P_mid/1e5:.4f} bar @ "
                      f"rho={rho_mid:.5f} kg/m^3")
            return rho_mid
        if P_mid < P_target_Pa: rho_min = rho_mid
        else:                   rho_max = rho_mid

    if verbose:
        print(f"Solver stopped (tol not reached): P={P_mid/1e5:.4f} bar "
              f"after 50 bisection steps.")
    return rho_mid


__all__ = ["invert_T_from_u", "find_rho_for_target_P"]
