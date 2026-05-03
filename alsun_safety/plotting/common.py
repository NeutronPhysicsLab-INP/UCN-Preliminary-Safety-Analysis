"""
Common plotting helpers for the AlSUN safety solver figures.

These reproduce the ``plot_sat_curves``, ``plot_isotherm_branches``, and
the Monte Carlo envelope plots used in the manuscript figures, but with
the per-figure styling parameters consolidated in :func:`set_publication_style`.
"""

from __future__ import annotations
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np


def set_publication_style(font_size: int = 12) -> None:
    """Apply manuscript-baseline plot styling globally."""
    plt.rcParams.update({
        "font.size":       font_size,
        "figure.autolayout": True,
        "axes.grid":       True,
        "grid.alpha":      0.3,
        "lines.linewidth": 2.0,
    })


def plot_sat_curve(db, ax=None):
    """Plot the saturation curve P_sat(T) of the loaded property database."""
    if ax is None:
        _, ax = plt.subplots(figsize=(7, 5), dpi=100)
    T = db.sat.T_grid
    P = db.sat.P_sat
    ax.plot(T, P / 1e5, "b-")
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel("Saturation pressure (bar)")
    ax.set_title("Deuterium saturation curve")
    return ax


def plot_isotherms(db, T_list: Iterable[float] = (30, 60, 100, 200, 300), ax=None):
    """Plot P(rho) for each loaded isotherm in ``T_list``."""
    if ax is None:
        _, ax = plt.subplots(figsize=(7, 5), dpi=100)
    for T in T_list:
        if T not in db.iso_at_T:
            continue
        atT = db.iso_at_T[T]
        for branch in (atT.liq, atT.gas):
            if branch is None:
                continue
            ax.loglog(branch.rho, branch.P / 1e5,
                      label=f"T = {T:.0f} K ({branch.phase})")
    ax.set_xlabel("Density (kg/m$^3$)")
    ax.set_ylabel("Pressure (bar)")
    ax.legend(loc="best", fontsize=9)
    ax.set_title("Deuterium isothermal P(rho)")
    return ax


def plot_mc_envelope(time, median, lo, hi, *, ax=None,
                     color: str = "C0", label: str = None,
                     band_label: str = "5th-95th percentile"):
    """Plot a Monte Carlo ensemble as median + percentile band."""
    if ax is None:
        _, ax = plt.subplots(figsize=(7, 5), dpi=100)
    line = ax.plot(time, median, color=color, label=label)[0]
    ax.fill_between(time, lo, hi, color=color, alpha=0.2,
                    label=band_label if label is not None else None)
    return ax


__all__ = ["set_publication_style", "plot_sat_curve", "plot_isotherms",
           "plot_mc_envelope"]
