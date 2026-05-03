"""Run a single deterministic LOVA transient WITHOUT the helium coolant gap.

Reproduces the counterfactual single-shell LOVA scenario (alternative
designs with a circulating deuterium loop and no helium buffer).

Run from the repository root:

    python scripts/run_lova_without_helium.py
"""

from __future__ import annotations

import matplotlib.pyplot as plt

from alsun_safety.config import default_config
from alsun_safety.eos import PropertyDB
from alsun_safety.scenarios.lova_no_helium import run_lova_without_helium_gap
from alsun_safety.plotting import set_publication_style


def main() -> None:
    cfg = default_config()
    set_publication_style()

    print("Loading NIST property database from data/...")
    db = PropertyDB.from_directory(
        sat_path=cfg.paths.saturation_file,
        isotherm_dir=cfg.paths.isotherm_dir,
    )

    print("Running LOVA transient without helium gap (~60 s wall clock)...")
    df = run_lova_without_helium_gap(db, cfg)

    fig, ax = plt.subplots(2, 1, figsize=(7, 7), dpi=100, sharex=True)
    ax[0].plot(df["t"], df["P1"] / 1e5, label="$P_{D2}$")
    ax[0].set_ylabel("Pressure (bar)")
    ax[0].legend()

    ax[1].plot(df["t"], df["T_wall"], label="$T_{wall}$")
    ax[1].plot(df["t"], df["T1"],     label="$T_{D2}$")
    ax[1].set_xlabel("Time (s)")
    ax[1].set_ylabel("Temperature (K)")
    ax[1].legend()

    fig.suptitle("LOVA without helium coolant gap (counterfactual)")
    fig.savefig("lova_without_helium.png", dpi=300)
    print("Saved: lova_without_helium.png")


if __name__ == "__main__":
    main()
