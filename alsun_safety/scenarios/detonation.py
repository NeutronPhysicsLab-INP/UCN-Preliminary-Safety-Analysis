"""
Detonation scenario — sequential layer-by-layer stress chain (Section 3.7).

Models a stoichiometric D2-air detonation as a bounding combustion case.
Pressure decays through three concentric volumes (LD2 chamber, He gap,
vacuum vessel) by isothermal ideal-gas expansion:

.. math::
    P_n = P_{n-1} \\frac{V_{n-1}}{V_{n-1} + V_n}

Each successive shell is then loaded by the corresponding pressure step
(P_in -> P_1 -> P_2). Stresses are compared against admissible limits
:math:`\\sigma_{adm} = \\sigma_{UTS} / \\eta` for both Aluminum 5056 and
the lead disk.
"""

from __future__ import annotations
import math
from dataclasses import dataclass, asdict
from typing import Dict, List

from ..config import SimConfig


def spherical_shell_segment_volume(R: float, t: float, theta_rad: float) -> float:
    """Volume of a spherical-shell segment (forward cap, thickness t)."""
    Ro = R + t
    return (2.0 * math.pi / 3.0) * (Ro ** 3 - R ** 3) * (1.0 - math.cos(theta_rad))


@dataclass
class DetonationResult:
    """One row of the detonation results table."""
    component: str
    load_step: str
    sigma_calc_MPa: float
    sigma_adm_MPa: float

    @property
    def status(self) -> str:
        return "SAFE" if self.sigma_calc_MPa < self.sigma_adm_MPa else "FAIL"

    def as_dict(self) -> Dict[str, str]:
        d = asdict(self)
        d["status"] = self.status
        return d


def run_detonation(cfg: SimConfig, R_curv: float,
                   *, V_vacuum_total: float = None) -> List[DetonationResult]:
    """Compute Von Mises stresses on the four shells under a bounding detonation.

    Parameters
    ----------
    cfg : :class:`SimConfig`
        Manuscript-baseline configuration.
    R_curv : float
        Spherical-cap radius (m), normally from :func:`solve_R_curv`.
    V_vacuum_total : float, optional
        Total free vacuum volume into which the detonation gas finally
        expands. If ``None``, computed from cfg.geometry.

    Returns
    -------
    List of :class:`DetonationResult` (one per shell).
    """
    geom = cfg.geometry
    deto = cfg.detonation
    al   = cfg.aluminum
    pb   = cfg.lead

    angle_rad = math.radians(geom.cap_angle_deg)

    # Admissible stresses
    sigma_adm_al = al.sigma_uts / al.eta_safety
    sigma_adm_pb = pb.sigma_uts / al.eta_safety  # lead uses same eta

    # Pressures
    P_in = deto.P_CJ * deto.DDT_factor

    # ------ Layer 1: inner spherical cap (LD2-bounding, thickness = inner-wall t) ------
    first_layer = P_in * R_curv / (2.0 * geom.t_inner_wall)

    # ------ Helium-gap expansion: P_1 ------
    V_inventory = geom.V_inventory
    V_he_cool = (
        spherical_shell_segment_volume(R_curv, geom.he_gap, angle_rad)
        + math.pi * geom.L_cyl * (
            (geom.R_outer_wall + geom.he_gap) ** 2 - geom.R_outer_wall ** 2
        )
    )
    P_1 = P_in * V_inventory / (V_inventory + V_he_cool)

    # ------ Layer 2: helium-jacket outer wall ------
    R_jacket_mid = geom.R_outer_wall + geom.he_gap + geom.t_he_jacket / 2.0
    second_layer = P_1 * R_jacket_mid / (2.0 * geom.t_he_jacket)

    # ------ Vacuum-vessel expansion: P_2 ------
    if V_vacuum_total is None:
        # Cylindrical-only estimate (matches notebook cell 63 form)
        V_vacuum_total = math.pi * (geom.R_vacuum_inner ** 2
                                    - (geom.R_outer_wall + geom.t_outer_wall) ** 2
                                    ) * geom.L_total

    P_2 = P_1 * (V_inventory + V_he_cool) / (V_inventory + V_he_cool + V_vacuum_total)

    # ------ Layer 3: vacuum vessel ------
    R_vac_mid = geom.R_outer_wall + geom.t_outer_wall + geom.t_vacuum / 2.0
    third_layer = P_2 * R_vac_mid / (2.0 * geom.t_vacuum)

    # ------ Layer 4: lead disk (plate bending under P_2) ------
    # Conservative bolt-spacing approximation
    a_lead = 0.384 * math.sqrt(2.0) / 2.0
    fourth_layer = 3.0 * P_2 * (a_lead ** 2) / (geom.t_lead ** 2)

    return [
        DetonationResult("Inner shell (Al)",  "P_in", first_layer / 1e6,  sigma_adm_al / 1e6),
        DetonationResult("Outer shell (Al)",  "P_1",  second_layer / 1e6, sigma_adm_al / 1e6),
        DetonationResult("Vacuum shell (Al)", "P_2",  third_layer / 1e6,  sigma_adm_al / 1e6),
        DetonationResult("Lead disk (Pb)",    "P_2",  fourth_layer / 1e6, sigma_adm_pb / 1e6),
    ]


__all__ = ["run_detonation", "DetonationResult", "spherical_shell_segment_volume"]
