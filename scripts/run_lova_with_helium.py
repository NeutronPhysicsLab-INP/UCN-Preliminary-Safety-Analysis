"""Run a single deterministic LOVA transient with the helium coolant gap.

Reproduces the deterministic baseline used to generate the LOVA figures
in §4.2 of the manuscript. For the Monte Carlo ensemble see the
:mod:`alsun_safety.uq` driver and the example pattern below.

Run from the repository root:

    python scripts/run_lova_with_helium.py
"""

from __future__ import annotations
import math

import matplotlib.pyplot as plt

from alsun_safety.config import default_config
from alsun_safety.eos import PropertyDB
from alsun_safety.scenarios.lova_helium import run_lova_with_helium_gap
from alsun_safety.plotting import set_publication_style


def _shell_masses(geom, rho_al: float = 2700.0):
    """Mass of inner-shell + outer (helium-jacket) shell, used for thermal inertia."""
    angle_rad = math.radians(geom.cap_angle_deg)
    factor = (1.0 - math.cos(angle_rad))

    # Inner shell (outer-D2 wall, t = t_outer_wall): cylindrical body + spherical cap
    R_in_outer = geom.R_outer_wall
    R_in_inner = R_in_outer - geom.t_outer_wall
    vol_in_cyl = math.pi * (R_in_outer ** 2 - R_in_inner ** 2) * geom.L_cyl
    vol_in_sph = (2.0 * math.pi / 3.0) * (R_in_outer ** 3 - R_in_inner ** 3) * factor
    mass_inner = (vol_in_cyl + vol_in_sph) * rho_al

    # Outer shell (helium-jacket outer wall, t = t_he_jacket)
    R_out_outer = geom.R_he_jacket + geom.t_he_jacket
    R_out_inner = geom.R_he_jacket
    vol_out_cyl = math.pi * (R_out_outer ** 2 - R_out_inner ** 2) * geom.L_cyl
    vol_out_sph = (2.0 * math.pi / 3.0) * (R_out_outer ** 3 - R_out_inner ** 3) * factor
    mass_outer = (vol_out_cyl + vol_out_sph) * rho_al

    return mass_inner, mass_outer


def main() -> None:
    cfg = default_config()
    set_publication_style()

    print("Loading NIST property database from data/...")
    db = PropertyDB.from_directory(
        sat_path=cfg.paths.saturation_file,
        isotherm_dir=cfg.paths.isotherm_dir,
    )

    mass_inner, mass_outer = _shell_masses(cfg.geometry, rho_al=cfg.aluminum.rho)
    print(f"Wall masses: inner = {mass_inner:.2f} kg,  outer = {mass_outer:.2f} kg")

    print("Running LOVA transient with helium gap (~60 s wall clock)...")
    df = run_lova_with_helium_gap(
        db, cfg,
        mass_inner=mass_inner, mass_outer=mass_outer,
    )

    fig, ax = plt.subplots(2, 1, figsize=(7, 7), dpi=100, sharex=True)
    ax[0].plot(df["t"], df["P1"] / 1e5,  label="$P_{D2}$")
    ax[0].plot(df["t"], df["P_he"] / 1e5, label="$P_{He}$")
    ax[0].set_ylabel("Pressure (bar)")
    ax[0].legend()

    ax[1].plot(df["t"], df["T_out"], label="$T_{wall, outer}$")
    ax[1].plot(df["t"], df["T_in"],  label="$T_{wall, inner}$")
    ax[1].plot(df["t"], df["T1"],    label="$T_{D2}$")
    ax[1].set_xlabel("Time (s)")
    ax[1].set_ylabel("Temperature (K)")
    ax[1].legend()

    fig.suptitle("LOVA with helium coolant gap")
    fig.savefig("lova_with_helium.png", dpi=300)
    print("Saved: lova_with_helium.png")


if __name__ == "__main__":
    main()
