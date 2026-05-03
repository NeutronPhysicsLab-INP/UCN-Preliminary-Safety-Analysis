"""Reproduce the R4-11 air-ingress calculation (manuscript Section 3.6).

Computes the free vacuum-vessel volume between the helium-jacket outer
wall and the vacuum-vessel inner wall, then converts to an asymptotic
air mass and an average mass-ingress rate over the 100 s pressurization
ramp.

Run from the repository root:

    python scripts/air_ingress.py
"""

from __future__ import annotations
import math

from alsun_safety.config import default_config
from alsun_safety.geometry import solve_R_curv, free_vacuum_volume


def main() -> None:
    cfg = default_config()
    geom = cfg.geometry
    deut = cfg.deuterium

    # 1. Solve for R_curv such that the LD2 inventory volume holds the design mass
    target_volume = deut.mass / deut.rho_liq_ref
    R_curv = solve_R_curv(geom, target_volume)
    print(f"R_curv (spherical-cap radius)  = {R_curv*1000:.1f} mm")

    # 2. Free vacuum volume
    fv = free_vacuum_volume(geom, R_curv)
    print(f"V_vac cylindrical annulus       = {fv.V_cyl_annulus*1e3:6.3f} L")
    print(f"V_vac spherical-cap shell       = {fv.V_cap_shell *1e3:6.3f} L")
    print(f"V_vac total                     = {fv.V_total     *1e3:6.3f} L")

    # 3. Mass-ingress rate
    P_atm = 101_325.0
    T_amb = cfg.lova.T_amb
    M_air = 28.96e-3
    R_univ = 8.31446
    rho_air = P_atm * M_air / (R_univ * T_amb)
    print(f"\nrho_air(300 K, 1 atm)           = {rho_air:.4f} kg/m^3")

    tau = cfg.lova.pressurization_t
    m_air_total = fv.V_total * rho_air
    m_dot_avg = m_air_total / tau
    print(f"Total air ingress               = {m_air_total*1e3:.3f} g")
    print(f"Average rate (tau = {tau:.0f} s)        = {m_dot_avg*1e6:.2f} mg/s "
          f"({m_dot_avg:.3e} kg/s)")

    # 4. Cross-check against the linear pressure ramp
    p_dot = (cfg.lova.p_vac_atm - cfg.lova.p_vac_initial) / cfg.lova.pressurization_t
    m_dot_check = p_dot * fv.V_total * M_air / (R_univ * T_amb)
    ratio = m_dot_check / m_dot_avg if m_dot_avg > 0 else float("nan")
    print(f"\nCross-check from p_dot~{p_dot:.0f} Pa/s :  "
          f"m_dot = {m_dot_check:.3e} kg/s  ({ratio:.3f}x of average)")


if __name__ == "__main__":
    main()
