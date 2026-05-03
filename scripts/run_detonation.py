"""Reproduce Table 3 of the manuscript: detonation structural demand vs. limits.

Run from the repository root:

    python scripts/run_detonation.py
"""

from __future__ import annotations

from alsun_safety.config import default_config
from alsun_safety.geometry import solve_R_curv
from alsun_safety.scenarios.detonation import run_detonation


def main() -> None:
    cfg = default_config()

    target_volume = cfg.deuterium.mass / cfg.deuterium.rho_liq_ref
    R_curv = solve_R_curv(cfg.geometry, target_volume)

    results = run_detonation(cfg, R_curv)

    print(f"{'Component':22s}{'Load':6s}  {'Calc.':>10s}  {'Limit':>10s}  Status")
    print("-" * 64)
    for r in results:
        print(f"{r.component:22s}{r.load_step:6s}  "
              f"{r.sigma_calc_MPa:>8.1f} MPa  "
              f"{r.sigma_adm_MPa:>8.1f} MPa  "
              f"{r.status}")


if __name__ == "__main__":
    main()
