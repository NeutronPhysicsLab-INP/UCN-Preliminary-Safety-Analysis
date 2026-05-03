"""
alsun_safety
============

Python solver for the analytical safety assessment of the liquid deuterium
premoderator container at the AlSUN ultra-cold neutron source (WWR-K
reactor). Released alongside the manuscript:

    *Preliminary Safety Assessment of Liquid Deuterium Premoderator
    Container for Ultra-Cold Neutron Source at WWR-K Reactor (AlSUN).*

Top-level subpackages
---------------------
- :mod:`alsun_safety.config`        Physical constants, geometry, scenario parameters
- :mod:`alsun_safety.geometry`      R_curv solver, free vacuum-volume calculation
- :mod:`alsun_safety.eos`           NIST-tabulated equation-of-state handler
- :mod:`alsun_safety.relief`        Rupture-disk + orifice flow models
- :mod:`alsun_safety.stress`        Thin-walled-vessel stress analysis
- :mod:`alsun_safety.scenarios`     Cooling-failure, LOVA, detonation simulators
- :mod:`alsun_safety.uq`            Monte Carlo uncertainty driver
- :mod:`alsun_safety.plotting`      Reusable plot helpers
"""

from .config import default_config, SimConfig
from .geometry import solve_R_curv, free_vacuum_volume

__version__ = "1.0.0"
__all__ = ["default_config", "SimConfig", "solve_R_curv", "free_vacuum_volume"]
