"""
Configuration / parameter sets for the AlSUN safety solver.

Every parameter the solver needs lives in a dataclass below. The default
values reproduce the configuration analysed in the accompanying manuscript
(Table 1 and Sections 3.5–3.7); override individual fields when scanning
parameter space.

Conventions
-----------
- SI units throughout unless the field name says otherwise.
- Lengths: m. Pressures: Pa. Temperatures: K. Energies: J. Power: W.
- Wall radii are *inner radii* of cylindrical shells (matching Table 1).
"""

from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------

R_UNIV = 8.314462618          # J / (mol K)
N_A    = 6.02214076e23        # 1 / mol
ATM_PA = 101_325.0            # Pa per atm
G      = 9.80665              # m / s^2


# ---------------------------------------------------------------------------
# Deuterium operational state
# ---------------------------------------------------------------------------

@dataclass
class Deuterium:
    """Reference operational state of the deuterium inventory."""
    M_molar:     float = 4.0282e-3      # kg/mol
    T_init:      float = 20.0           # K, normal-operation temperature
    p_init:      float = 1.5 * ATM_PA   # Pa, normal-operation pressure
    mass:        float = 14.0           # kg, design LD2 inventory
    rho_liq_ref: float = 169.0          # kg/m^3, liquid density near 20 K
    rho_gas_ref: float = 2.105          # kg/m^3, vapor density at saturation
    rho_gas_300: float = 0.1645         # kg/m^3, deuterium gas at 300 K, 1 atm
    T_triple:    float = 18.7           # K, triple point
    T_boil:      float = 25.1           # K, normal boiling point
    cp_liq_ref:  float = 6.3e3          # J/(kg K), reference liquid cp
    cp_gas_ref:  float = 5.1847e3       # J/(kg K), reference gas cp


# ---------------------------------------------------------------------------
# Geometry — Table 1 of the manuscript
# ---------------------------------------------------------------------------

@dataclass
class Geometry:
    """Cylindrical-with-spherical-cap multi-shell deuterium chamber.

    All lengths in metres. Wall radii are inner radii (matching Table 1).
    Wall thicknesses are wall-only (not gap distances).
    """
    # ------ Deuterium chamber (LD2-bounding shells) ------
    R_inner_wall:   float = 0.160       # inner deuterium wall, R = 160 mm
    t_inner_wall:   float = 0.003       # 3 mm Al 5056

    R_outer_wall:   float = 0.230       # outer deuterium wall, R = 230 mm
    t_outer_wall:   float = 0.003       # 3 mm Al 5056

    # ------ Helium cooling envelope ------
    he_gap:         float = 0.002       # 2 mm He coolant gap (between R_outer_wall outer face
                                         # and helium-jacket outer wall inner face)
    R_he_jacket:    float = 0.235       # helium-jacket outer wall, R = 235 mm
    t_he_jacket:    float = 0.004       # 4 mm Al 5056

    # ------ Vacuum vessel ------
    R_vacuum_inner: float = 0.243       # vacuum-vessel inner radius, R = 243 mm
    t_vacuum:       float = 0.009       # 9 mm Al 5056

    # ------ Lead shielding disk + annular plate (rear) ------
    R_lead:         float = 0.485       # lead disk radius, R = 485 mm
    t_lead:         float = 0.100       # 100 mm Pb

    R_plate:        float = 0.230       # annular plate radius (rear, sealing the LD2 chamber)
    t_plate:        float = 0.015       # 15 mm Al 5056 (omitted from present stress evaluation)

    # ------ Cylindrical body length + spherical cap geometry ------
    L_cyl:          float = 0.510       # cylindrical body length, m
    L_total:        float = 0.710       # total assembly length (incl. cap + plate)
    AC_intercept:   float = 0.200       # m — spherical-cap geometry intercept (cell 1 of notebook)
    cap_angle_deg:  float = 80.0        # spherical-cap half-angle (degrees)

    # ------ Inventory volume (set so 14 kg LD2 fits at rho_liq) ------
    V_inventory:    float = 0.085       # m^3 = 85 L (matches manuscript / notebook)


# ---------------------------------------------------------------------------
# Materials
# ---------------------------------------------------------------------------

@dataclass
class Aluminum5056:
    """Aluminum Alloy 5056 — wall material for all aluminium shells."""
    name:           str   = "Al 5056"
    rho:            float = 2700.0      # kg/m^3
    k_thermal:      float = 160.0       # W/(m K)
    sigma_yield_02: float = 150e6       # Pa, 0.2 % proof stress (cited)
    sigma_uts:      float = 290e6       # Pa, ultimate tensile strength (cited)
    eta_safety:     float = 1.5         # Safety factor (ASME BPVC §VIII Div.1)


@dataclass
class Helium:
    """Helium coolant in the helium-jacket envelope."""
    p_init:    float = 1.0 * ATM_PA     # Pa
    T_init:    float = 20.0             # K
    R_specific: float = 2077.0          # J/(kg K), specific gas constant
    k_thermal: float = 0.0267           # W/(m K) at 20 K, 1 atm


@dataclass
class Lead:
    """Lead shielding disk."""
    name:      str   = "Pb"
    rho:       float = 11340.0          # kg/m^3
    sigma_uts: float = 15e6             # Pa, room-temperature UTS (ASM Handbook)


# ---------------------------------------------------------------------------
# Rupture-disk + ballast-tank relief
# ---------------------------------------------------------------------------

@dataclass
class RuptureDisk:
    """OsecoElfab LoKr-class rupture disk + warm ballast tank."""
    p_burst:        float = 2.0 * ATM_PA  # Pa, absolute setpoint
    burst_tol:      float = 0.05          # ±5 %, manufacturer spec (uniform)
    diameter:       float = 0.025         # m (25 mm)
    Cd:             float = 0.8           # Discharge coefficient (conservative)
    K_R:            float = 0.22          # Flow resistance factor (manufacturer)
    open_time:      float = 0.1           # s — characteristic disk-opening time
    mdot_max_cap:   float = 50.0          # kg/s — pipe choke ceiling

    # Ballast-tank inventory (warm receiving volume for blowdown)
    V_ballast:      float = 5.0           # m^3
    T_ballast_init: float = 300.0         # K
    p_ballast_init: float = 1.0 * ATM_PA  # Pa


# ---------------------------------------------------------------------------
# Cooling-failure scenario heat budget (Section 3.5)
# ---------------------------------------------------------------------------

@dataclass
class CoolingFailureLoads:
    """Total nominal heat load to the deuterium inventory under loss of helium cooling.

    Values match the manuscript Section 3.5 (24 + 69.1 + 78.4 = 171.5 W).
    The nuclear-heating contributions are MCNP outputs for the AlSUN baseline
    configuration with Aluminum Alloy 5056 vessels; they are reported in the
    accompanying paper and are reproduced here verbatim so the published
    figures can be regenerated from this code.
    """
    Q_radiative:   float = 24.0           # W — radiative input (geometry-specific)
    Q_nuclear_D2:  float = 69.1           # W — MCNP, deposition in deuterium inventory
    Q_nuclear_Al:  float = 78.4           # W — MCNP, deposition in aluminium structure

    @property
    def Q_total(self) -> float:
        return self.Q_radiative + self.Q_nuclear_D2 + self.Q_nuclear_Al


# ---------------------------------------------------------------------------
# LOVA scenario (Section 3.6)
# ---------------------------------------------------------------------------

@dataclass
class LOVA:
    """Loss-of-Vacuum Accident: air ingress + cryopumping deposition."""
    p_vac_initial:    float = 1.3e-3      # Pa — initial high vacuum
    p_vac_atm:        float = ATM_PA      # Pa — final atmospheric
    pressurization_t: float = 100.0       # s — linear ramp duration
    flux_deposition:  float = 38e3        # W/m^2 — bare-vessel cryopumping flux (Lehmann/Zahn)
    P_half_knudsen:   float = 100.0       # Pa — Knudsen-flow half-conductivity (Corruccini)
    h_natural_conv:   float = 15.0        # W/(m^2 K) — Churchill-Chu order of magnitude
    h_nucleate_boil:  float = 12000.0     # W/(m^2 K) — Drayer & Timmerhaus
    h_film_boil:      float = 500.0       # W/(m^2 K) — Wang et al.
    epsilon_radiative: float = 0.055      # Polished Al cryogenic emissivity (Barron & Nellis)
    T_amb:            float = 298.0       # K — ambient


# ---------------------------------------------------------------------------
# Detonation scenario (Section 3.7)
# ---------------------------------------------------------------------------

@dataclass
class Detonation:
    """Stoichiometric H2/D2-air detonation as bounding combustion case."""
    P_CJ_bar:     float = 15.6            # bar — Chapman–Jouguet pressure (Shepherd 2009)
    P_CJ:         float = 15.6e5          # Pa
    DDT_factor:   float = 5.0             # Reflected-shock + DDT amplification, fully confined
    U_CJ:         float = 1971.0          # m/s — CJ velocity (Shepherd 2009)
    lambda_cell:  float = 0.30            # m — critical tube diameter (Koroll 1991)


# ---------------------------------------------------------------------------
# Monte Carlo / UQ
# ---------------------------------------------------------------------------

@dataclass
class UQ:
    """Uncertainty-quantification driver settings."""
    n_runs:         int   = 100
    sim_time:       float = 30000.0       # s — cooling-failure ensemble length
    sim_time_lova:  float = 5000.0        # s — LOVA ensemble length
    decimation:     int   = 10            # save every Nth row
    n_processes:    int   = 12            # multiprocessing.Pool size

    # Sampled-input distributions (matches Section 3.8)
    burst_rel_halfwidth: float = 0.05     # ±5 % uniform on rupture-disk setpoint
    Q_rel_sigma:         float = 0.10     # σ = 10 % Gaussian on cooling-failure heat load
    mass_rel_sigma:      float = 0.005    # σ = 0.5 % Gaussian on initial D2 mass

    seed:           Optional[int] = None  # reseeded per worker process


# ---------------------------------------------------------------------------
# Property-data paths
# ---------------------------------------------------------------------------

@dataclass
class PropertyPaths:
    """Paths to the NIST tabular property files (place under data/)."""
    data_dir:         Path = field(default_factory=lambda: Path("data"))
    saturation_file:  Path = field(default_factory=lambda: Path("data") / "fluid-2.cgi")
    isotherm_dir:     Path = field(default_factory=lambda: Path("data") / "isotherms")


# ---------------------------------------------------------------------------
# Solver / numerics
# ---------------------------------------------------------------------------

@dataclass
class SolverOptions:
    """Numerical solver options (time-step adaptation, tolerances)."""
    dt_initial:      float = 0.01         # s
    dt_burst:        float = 0.0005       # s — fine step around disk burst
    dt_quasi:        float = 0.1          # s — coarse step in quasi-steady regime
    eos_tol:         float = 1e-3         # J/kg — energy-inversion tolerance
    eos_max_iter:    int   = 200
    output_csv:      str   = "results.csv"


# ---------------------------------------------------------------------------
# Top-level configuration
# ---------------------------------------------------------------------------

@dataclass
class SimConfig:
    """Bundle of every configuration block.

    Construct via :func:`default_config` for the manuscript baseline; override
    individual fields when scanning parameter space.
    """
    deuterium:         Deuterium             = field(default_factory=Deuterium)
    geometry:          Geometry              = field(default_factory=Geometry)
    aluminum:          Aluminum5056          = field(default_factory=Aluminum5056)
    helium:            Helium                = field(default_factory=Helium)
    lead:              Lead                  = field(default_factory=Lead)
    rupture_disk:      RuptureDisk           = field(default_factory=RuptureDisk)
    cooling_failure:   CoolingFailureLoads   = field(default_factory=CoolingFailureLoads)
    lova:              LOVA                  = field(default_factory=LOVA)
    detonation:        Detonation            = field(default_factory=Detonation)
    uq:                UQ                    = field(default_factory=UQ)
    paths:             PropertyPaths         = field(default_factory=PropertyPaths)
    solver:            SolverOptions         = field(default_factory=SolverOptions)


def default_config() -> SimConfig:
    """Return the manuscript-baseline configuration."""
    return SimConfig()


__all__ = [
    "ATM_PA", "G", "N_A", "R_UNIV",
    "Deuterium", "Geometry", "Aluminum5056", "Helium", "Lead",
    "RuptureDisk", "CoolingFailureLoads", "LOVA", "Detonation",
    "UQ", "PropertyPaths", "SolverOptions",
    "SimConfig", "default_config",
]
