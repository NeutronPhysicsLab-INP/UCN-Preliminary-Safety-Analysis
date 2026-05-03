"""
Thin-walled-vessel stress analysis.

Implements the analytical models of §3.3 of the manuscript:

- :func:`calculate_dynamic_stresses` — single-shell vessel: cylindrical
  body, spherical forward cap, and rear annular plate. Used for the
  cooling-failure scenario where there is no helium-gap pressurization.
- :func:`calculate_2wall_stresses` — two-wall (deuterium chamber + helium
  jacket) Von Mises stresses, used for the LOVA scenario where the
  trapped helium pressurizes both shells differentially.

Both routines accept a pandas DataFrame of transient pressures
(``df["P1"]`` in Pa, optionally ``df["P_he"]``) and return a DataFrame of
Von Mises stresses in MPa.
"""

from __future__ import annotations
import numpy as np
import pandas as pd


def calculate_dynamic_stresses(df: pd.DataFrame,
                               R_cyl:    float, t_cyl:    float,
                               R_sphere: float, t_sphere: float,
                               R_plate:  float, t_plate:  float,
                               kappa: float = 0.273) -> pd.DataFrame:
    """Compute Von Mises stress histories for cylinder, sphere, and plate.

    The cylindrical and spherical components use the Mariotte (thin-walled)
    membrane stresses with an averaged radial component; the rear annular
    plate is treated as a simply-supported circular plate under uniform
    pressure with bending coefficient :math:`\\kappa`.

    Parameters
    ----------
    df : DataFrame
        Must contain ``df["t"]`` (s) and ``df["P1"]`` (Pa).
    R_cyl, t_cyl : cylindrical body radius and wall thickness (m).
    R_sphere, t_sphere : forward cap radius and thickness (m).
    R_plate, t_plate : rear plate radius and thickness (m).
    kappa : plate bending coefficient (default 0.273, simply-supported
        circular plate, central deflection per Timoshenko & Woinowsky-Krieger).

    Returns
    -------
    DataFrame with columns ``t``, ``Stress_Cyl_MPa``, ``Stress_Sphere_MPa``,
    ``Stress_Plate_MPa``.
    """
    P = df["P1"].to_numpy()  # Pa

    # 1. Cylindrical wall: hoop, axial, radial -> Von Mises
    s_theta = P * R_cyl / t_cyl
    s_z     = P * R_cyl / (2.0 * t_cyl)
    s_r     = -P / 2.0
    vm_cyl = np.sqrt(0.5 * (
        (s_theta - s_z) ** 2
        + (s_theta - s_r) ** 2
        + (s_r - s_z) ** 2
    ))

    # 2. Spherical cap: sigma_theta = sigma_phi
    s_sphere   = P * R_sphere / (2.0 * t_sphere)
    s_r_sphere = -P / 2.0
    vm_sphere  = np.abs(s_sphere - s_r_sphere)  # reduces to |s_sphere - s_r|

    # 3. Annular plate: bending stress
    vm_plate = kappa * P * (R_plate / t_plate) ** 2

    return pd.DataFrame({
        "t": df["t"].to_numpy(),
        "Stress_Cyl_MPa":    vm_cyl    / 1e6,
        "Stress_Sphere_MPa": vm_sphere / 1e6,
        "Stress_Plate_MPa":  vm_plate  / 1e6,
    })


def calculate_2wall_stresses(df: pd.DataFrame,
                             R_in:  float, t_in:  float,
                             R_out: float, t_out: float,
                             P_atm: float = 101_325.0) -> pd.DataFrame:
    """Two-wall LOVA stress analysis (outer-D2 wall + helium-jacket outer wall).

    During a LOVA the helium-gap pressure rises above both the deuterium
    pressure and atmospheric, putting the outer-D2 wall in compression and
    the helium-jacket outer wall in tension. The Von Mises stresses on both
    walls are computed from the differential pressure across each wall and
    returned with sign attached (sign(dP)) so that buckling vs. yielding
    can be distinguished downstream.

    Parameters
    ----------
    df : DataFrame
        Must contain ``df["t"]``, ``df["P1"]`` (LD2 pressure, Pa) and
        ``df["P_he"]`` (helium-gap pressure, Pa).
    R_in, t_in   : inner wall (LD2-He boundary) radius and thickness (m).
    R_out, t_out : outer wall (He-vacuum boundary) radius and thickness (m).
    P_atm : ambient pressure (Pa).
    """
    P_d2 = df["P1"].to_numpy()
    P_he = df["P_he"].to_numpy()

    # Inner wall sees (P_he - P_D2): tensile if D2 > He, compressive if He > D2
    dP_inner = P_d2 - P_he
    s_h_in = dP_inner * R_in / t_in
    s_z_in = dP_inner * R_in / (2.0 * t_in)
    s_r_in = -np.abs(dP_inner) / 2.0
    vm_in = np.sqrt(0.5 * (
        (s_h_in - s_z_in) ** 2
        + (s_h_in - s_r_in) ** 2
        + (s_r_in - s_z_in) ** 2
    ))

    # Outer wall sees (P_he - P_atm), always tensile during a LOVA
    dP_outer = P_he - P_atm
    s_h_out = dP_outer * R_out / t_out
    s_z_out = dP_outer * R_out / (2.0 * t_out)
    s_r_out = -dP_outer / 2.0
    vm_out = np.sqrt(0.5 * (
        (s_h_out - s_z_out) ** 2
        + (s_h_out - s_r_out) ** 2
        + (s_r_out - s_z_out) ** 2
    ))

    return pd.DataFrame({
        "t": df["t"].to_numpy(),
        "P_D2_bar":          P_d2 / 1e5,
        "P_He_bar":          P_he / 1e5,
        # Sign: + tensile, - compressive
        "Stress_Inner_MPa":  vm_in  / 1e6 * np.sign(dP_inner),
        "Stress_Outer_MPa":  vm_out / 1e6,
    })


__all__ = ["calculate_dynamic_stresses", "calculate_2wall_stresses"]
