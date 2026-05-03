"""Thin-walled-vessel stress analysis."""
from .thin_wall import calculate_dynamic_stresses, calculate_2wall_stresses

__all__ = ["calculate_dynamic_stresses", "calculate_2wall_stresses"]
