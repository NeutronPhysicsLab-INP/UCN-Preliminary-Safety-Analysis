"""
Geometry computations for the AlSUN deuterium chamber.

Two utilities used by the scenarios:

- :func:`solve_R_curv` — bisection solver for the spherical-cap radius
  ``R_curv`` such that the total interior volume of the cylindrical-with-
  spherical-cap chamber equals the LD2 inventory volume at the design
  density. (Replaces the ad-hoc bisection loop in cell 1 of the original
  notebook.)

- :func:`free_vacuum_volume` — analytical estimate of the air-ingress
  volume between the helium-jacket outer wall and the vacuum-vessel inner
  wall. Used by the LOVA mass-ingress calculation in §3.6 of the
  manuscript and in :file:`scripts/air_ingress.py`.
"""

from __future__ import annotations
import math
from dataclasses import dataclass
from typing import Tuple

from .config import Geometry


def total_volume(R_curv: float, R_out: float, R_in: float,
                 L_cyl: float, AC: float, cap_angle_deg: float) -> Tuple[float, float]:
    """Total LD2 chamber interior volume + cap-only volume.

    Returns
    -------
    (V_total, V_cap) in m^3.
    """
    angle = cap_angle_deg
    OB = math.sqrt(R_curv ** 2 - R_out ** 2)
    AB = R_curv - OB
    BC = AC - AB
    V_cap = (
        (4.0 / 3.0) * math.pi * R_curv ** 3 * (angle / 360.0)
        - (1.0 / 3.0) * math.pi * (R_out ** 2) * OB
    )
    V_total = (
        V_cap
        + math.pi * (R_out ** 2) * BC
        + math.pi * (R_out ** 2 - R_in ** 2) * L_cyl
    )
    return V_total, V_cap


def solve_R_curv(geom: Geometry, target_volume: float,
                 R_lower: float = 0.20, R_upper: float = 1.00,
                 tol: float = 1e-2, max_iter: int = 1000) -> float:
    """Find R_curv such that the chamber volume equals ``target_volume``.

    Bisection, with the same tolerance used in cell 1 of the original
    notebook (1 % of target volume).
    """
    angle = geom.cap_angle_deg
    R_out = geom.R_outer_wall
    R_in  = geom.R_inner_wall
    L_cyl = geom.L_cyl
    AC    = geom.AC_intercept

    for _ in range(max_iter):
        R_mid = 0.5 * (R_lower + R_upper)
        V, _  = total_volume(R_mid, R_out, R_in, L_cyl, AC, angle)
        if abs(V - target_volume) < tol * target_volume:
            return math.ceil(R_mid * 100) / 100  # round up to nearest cm
        if V < target_volume:
            R_lower = R_mid
        else:
            R_upper = R_mid
    # Did not converge: return the best midpoint anyway
    return math.ceil(0.5 * (R_lower + R_upper) * 100) / 100


@dataclass
class FreeVacuumVolume:
    """Result bundle from :func:`free_vacuum_volume`."""
    V_cyl_annulus: float    # cylindrical annulus contribution (m^3)
    V_cap_shell:   float    # spherical-cap shell contribution (m^3)
    V_total:       float    # sum (m^3)


def free_vacuum_volume(geom: Geometry, R_curv: float,
                       gap_sph: float = 0.015,
                       cap_wall_to_vacuum_offset: float = None,
                       ) -> FreeVacuumVolume:
    """Free vacuum-vessel volume between helium-jacket outer face and
    vacuum-vessel inner face.

    The geometry is asymmetric between the cylindrical body and the
    spherical front cap:

    - **Cylindrical annulus**: the radial gap is fixed by the published
      Table 1 radii, ``R_vacuum_inner - (R_he_jacket + t_he_jacket)``.

    - **Spherical-cap shell**: the inner shell surface is the helium-jacket
      cap outer face at ``R_curv + cap_wall_to_vacuum_offset`` (sum of the
      cap-region wall thicknesses + the He gap; default 7 mm =
      3 mm inner cap + 2 mm He + 2 mm outer cap). The outer shell surface
      is the vacuum-vessel cap inner face at
      ``R_curv + cap_wall_to_vacuum_offset + gap_sph`` (default 15 mm
      vacuum gap in the cap region; this is the value confirmed by the
      AlSUN design).

    Parameters
    ----------
    geom    : :class:`Geometry` instance with the as-built dimensions.
    R_curv  : LD2-cap spherical radius (m), normally from :func:`solve_R_curv`.
    gap_sph : radial vacuum gap in the spherical-cap region (m). Default
              0.015 m (15 mm), per the manuscript description and cell 32
              of the source notebook.
    cap_wall_to_vacuum_offset : distance from the LD2 cap surface to the
              helium-jacket cap outer face. Default 0.007 m, taken from
              the cap-region wall thicknesses (3 + 2 + 2 mm). Override if
              your cap-region wall thicknesses differ from the spherical
              defaults.
    """
    R_he_jacket_outer_cyl = geom.R_he_jacket + geom.t_he_jacket
    R_vac_inner_cyl       = geom.R_vacuum_inner
    L_cyl                 = geom.L_cyl

    # Cylindrical annulus (Table 1 radii)
    V_cyl = math.pi * (R_vac_inner_cyl ** 2 - R_he_jacket_outer_cyl ** 2) * L_cyl

    # Cap-region inner / outer surfaces of the vacuum shell.
    # In the cylindrical body the outer wall is 4 mm (cfg.geometry.t_he_jacket)
    # but on the front cap it is 2 mm (cell 32 of the original notebook,
    # parameter t_out_sph). Hence the cap-region offset uses the cap thicknesses,
    # not the cylindrical ones:
    #     3 mm (outer-D2 cap wall, t_in_sph in cell 32)
    #   + 2 mm (He gap, t_he)
    #   + 2 mm (helium-jacket cap outer wall, t_out_sph in cell 32)
    #   = 7 mm
    if cap_wall_to_vacuum_offset is None:
        cap_wall_to_vacuum_offset = 0.007

    R_inner_cap = R_curv + cap_wall_to_vacuum_offset       # he-jacket cap outer face
    R_outer_cap = R_inner_cap + gap_sph                    # vacuum-vessel cap inner face

    angle_rad = math.radians(geom.cap_angle_deg)
    V_cap = (2.0 * math.pi / 3.0) * (R_outer_cap ** 3 - R_inner_cap ** 3) \
            * (1.0 - math.cos(angle_rad))

    return FreeVacuumVolume(
        V_cyl_annulus=V_cyl,
        V_cap_shell=V_cap,
        V_total=V_cyl + V_cap,
    )


__all__ = ["total_volume", "solve_R_curv", "free_vacuum_volume", "FreeVacuumVolume"]
