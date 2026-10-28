"""Effective diffusion coefficient of hydrophilic molecules in lipidic cubic phases."""

from .core import (CUBIC_GEOMETRY, PHASES, DiffusionResult, cubic_lipid_length, diffusion,
                   reduced_diffusion, rescale_D0, water_viscosity, water_volume_fraction)

__version__ = "1.0.0"

__all__ = ["CUBIC_GEOMETRY", "PHASES", "DiffusionResult", "cubic_lipid_length", "diffusion",
           "reduced_diffusion", "rescale_D0", "water_viscosity", "water_volume_fraction"]
