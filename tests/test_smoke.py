"""
Lightweight smoke tests that exercise the no-NIST-data code paths.

These run without the property database and confirm that the package
imports cleanly, the geometry solver converges, V_vac matches the
manuscript, and the detonation analysis produces the expected pass/fail
pattern. Run with::

    pytest tests/

or simply::

    python tests/test_smoke.py
"""

from __future__ import annotations

import math


def test_package_imports():
    """Every submodule should import cleanly."""
    import alsun_safety
    from alsun_safety import config, geometry
    from alsun_safety.eos       import property_db, solver
    from alsun_safety.relief    import flow, rupture_disk
    from alsun_safety.stress    import thin_wall
    from alsun_safety.scenarios import (
        cooling_failure, lova_helium, lova_no_helium, detonation,
    )
    from alsun_safety.uq        import monte_carlo
    from alsun_safety.plotting  import common
    assert alsun_safety.__version__ == "1.0.0"


def test_R_curv_solver():
    """R_curv from inventory-volume constraint must be 340 mm."""
    from alsun_safety.config import default_config
    from alsun_safety.geometry import solve_R_curv
    cfg = default_config()
    target_volume = cfg.deuterium.mass / cfg.deuterium.rho_liq_ref
    R_curv = solve_R_curv(cfg.geometry, target_volume)
    assert math.isclose(R_curv, 0.340, abs_tol=0.01), \
        f"R_curv expected 0.340 m, got {R_curv:.4f}"


def test_free_vacuum_volume():
    """V_vac must reproduce the manuscript value (≈12.9 L)."""
    from alsun_safety.config import default_config
    from alsun_safety.geometry import solve_R_curv, free_vacuum_volume
    cfg = default_config()
    target_volume = cfg.deuterium.mass / cfg.deuterium.rho_liq_ref
    R_curv = solve_R_curv(cfg.geometry, target_volume)
    fv = free_vacuum_volume(cfg.geometry, R_curv)

    assert math.isclose(fv.V_cyl_annulus * 1e3, 3.09, abs_tol=0.05), \
        f"Expected V_cyl_annulus ≈ 3.09 L, got {fv.V_cyl_annulus*1e3:.3f} L"
    assert math.isclose(fv.V_cap_shell * 1e3, 9.79, abs_tol=0.10), \
        f"Expected V_cap_shell ≈ 9.79 L, got {fv.V_cap_shell*1e3:.3f} L"
    assert math.isclose(fv.V_total * 1e3, 12.88, abs_tol=0.15), \
        f"Expected V_total ≈ 12.88 L, got {fv.V_total*1e3:.3f} L"


def test_detonation_results():
    """Detonation-analysis must reproduce the manuscript Table 3 pass/fail
    pattern: inner & outer Al shells FAIL, vacuum shell SAFE, lead disk FAIL."""
    from alsun_safety.config import default_config
    from alsun_safety.geometry import solve_R_curv
    from alsun_safety.scenarios.detonation import run_detonation
    cfg = default_config()
    target_volume = cfg.deuterium.mass / cfg.deuterium.rho_liq_ref
    R_curv = solve_R_curv(cfg.geometry, target_volume)
    results = run_detonation(cfg, R_curv)

    statuses = {r.component: r.status for r in results}
    assert statuses["Inner shell (Al)"]  == "FAIL"
    assert statuses["Outer shell (Al)"]  == "FAIL"
    assert statuses["Vacuum shell (Al)"] == "SAFE"
    assert statuses["Lead disk (Pb)"]    == "FAIL"


def test_air_ingress_rate():
    """Air-ingress rate ≈ 152 mg/s for the manuscript geometry."""
    from alsun_safety.config import default_config
    from alsun_safety.geometry import solve_R_curv, free_vacuum_volume
    cfg = default_config()
    target_volume = cfg.deuterium.mass / cfg.deuterium.rho_liq_ref
    R_curv = solve_R_curv(cfg.geometry, target_volume)
    fv = free_vacuum_volume(cfg.geometry, R_curv)

    P_atm   = 101325.0
    T_amb   = cfg.lova.T_amb
    M_air   = 28.96e-3
    R_univ  = 8.31446
    rho_air = P_atm * M_air / (R_univ * T_amb)
    m_dot   = fv.V_total * rho_air / cfg.lova.pressurization_t

    assert math.isclose(m_dot * 1e6, 152.5, abs_tol=2.0), \
        f"Expected m_dot ≈ 152 mg/s, got {m_dot*1e6:.2f} mg/s"


if __name__ == "__main__":
    # Run all tests when invoked directly
    test_package_imports()
    test_R_curv_solver()
    test_free_vacuum_volume()
    test_detonation_results()
    test_air_ingress_rate()
    print("All smoke tests passed.")
