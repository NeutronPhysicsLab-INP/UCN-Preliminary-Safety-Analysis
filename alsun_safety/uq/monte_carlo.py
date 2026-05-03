"""
Generic Monte Carlo driver for uncertainty quantification.

Re-implements the multiprocessing pattern used in the cooling-failure and
LOVA UQ blocks of the original notebook in a scenario-agnostic form. The
caller provides a ``run_one(sample) -> DataFrame`` callable; this module
generates ``N`` perturbed input samples, runs them in parallel, and
returns the concatenated DataFrame.

The default perturbation distributions match Section 3.8 of the
manuscript: ±5 % uniform on the rupture-disk burst pressure, σ = 10 %
Gaussian on the heat load, σ = 0.5 % Gaussian on the deuterium mass.
"""

from __future__ import annotations
import multiprocessing
import time
from dataclasses import dataclass
from typing import Callable, Iterable, List, Tuple

import numpy as np
import pandas as pd


@dataclass
class MCSample:
    """One Monte Carlo input perturbation."""
    run_id:     int
    m_factor:   float    # multiplier on nominal LD2 mass
    P_factor:   float    # multiplier on nominal burst pressure
    Q_factor:   float    # multiplier on nominal heat load


def draw_samples(n_runs: int,
                 mass_rel_sigma:      float = 0.005,
                 burst_rel_halfwidth: float = 0.05,
                 Q_rel_sigma:         float = 0.10,
                 *, seed: int = None) -> List[MCSample]:
    """Generate ``n_runs`` perturbation samples (Section 3.8 distributions)."""
    rng = np.random.default_rng(seed)
    return [
        MCSample(
            run_id   = i + 1,
            m_factor = float(rng.normal(1.0, mass_rel_sigma)),
            P_factor = float(rng.uniform(1.0 - burst_rel_halfwidth,
                                         1.0 + burst_rel_halfwidth)),
            Q_factor = float(rng.normal(1.0, Q_rel_sigma)),
        )
        for i in range(n_runs)
    ]


def run_ensemble(run_one: Callable[[MCSample], pd.DataFrame],
                 samples: Iterable[MCSample],
                 *,
                 n_processes: int = 12,
                 verbose: bool = True) -> List[pd.DataFrame]:
    """Run ``run_one`` on each sample in parallel; return list of per-run DataFrames.

    The ``run_one`` callable must be importable at module scope (a quirk
    of multiprocessing.Pool). Each returned DataFrame should already be
    decimated to manageable size; the driver does not re-decimate.
    """
    samples = list(samples)
    t0 = time.time()
    if verbose:
        print(f"[MC] Starting {len(samples)} runs on {n_processes} processes...")

    with multiprocessing.Pool(processes=n_processes) as pool:
        results = pool.map(run_one, samples)

    if verbose:
        print(f"[MC] Finished in {time.time() - t0:.1f} s.")
    return results


def aggregate(results: List[pd.DataFrame],
              *, common_time: np.ndarray,
              key: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Interpolate ``key`` from each DataFrame onto ``common_time`` and stack.

    Returns
    -------
    median, p5, p95 : ndarray, each shaped like ``common_time``.
    """
    from scipy.interpolate import interp1d

    rows = []
    for df in results:
        t = df["t"].to_numpy()
        y = df[key].to_numpy()
        f = interp1d(t, y, bounds_error=False, fill_value=(y[0], y[-1]))
        rows.append(f(common_time))
    M = np.vstack(rows)
    return (np.median(M, axis=0),
            np.percentile(M, 5, axis=0),
            np.percentile(M, 95, axis=0))


__all__ = ["MCSample", "draw_samples", "run_ensemble", "aggregate"]
