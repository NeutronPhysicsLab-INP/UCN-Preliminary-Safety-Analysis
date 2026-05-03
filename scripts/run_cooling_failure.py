"""Run a single deterministic cooling-failure transient and plot the result.

For the ensemble (Monte Carlo) version see ``scripts/run_cooling_failure_uq.py``
which is left as a small wrapper around this entry point.

Run from the repository root:

    python scripts/run_cooling_failure.py
"""

from __future__ import annotations
import matplotlib.pyplot as plt

from alsun_safety.config import default_config
from alsun_safety.eos import PropertyDB
from alsun_safety.scenarios.cooling_failure import run_cooling_failure
from alsun_safety.plotting import set_publication_style


def main() -> None:
    cfg = default_config()
    set_publication_style()

    print("Loading NIST property database from data/...")
    db = PropertyDB.from_directory(
        sat_path=cfg.paths.saturation_file,
        isotherm_dir=cfg.paths.isotherm_dir,
    )

    print("Running cooling-failure transient (~30 s wall clock)...")
    hist = run_cooling_failure(db, cfg)

    fig, ax = plt.subplots(2, 1, figsize=(7, 7), dpi=100, sharex=True)
    ax[0].plot(hist["t"], hist["P1"], label="Vessel ($P_1$)")
    ax[0].plot(hist["t"], hist["P2"], label="Ballast ($P_2$)")
    ax[0].set_ylabel("Pressure (bar)")
    ax[0].legend()

    ax[1].plot(hist["t"], hist["m1"], label="Vessel ($m_1$)")
    ax[1].plot(hist["t"], hist["m2"], label="Ballast ($m_2$)")
    ax[1].set_xlabel("Time (s)")
    ax[1].set_ylabel("Mass (kg)")
    ax[1].legend()

    fig.suptitle("Cooling-failure transient")
    fig.savefig("cooling_failure.png", dpi=300)
    print("Saved: cooling_failure.png")


if __name__ == "__main__":
    main()
