"""Effective diffusion coefficient of hydrophilic molecules in lipidic cubic phases.

Python port of the "Diffusion Coefficient in Water Channels" tool available at
https://www.salvatoreassenza.com/#tools, following the analysis of
S. Assenza and R. Mezzenga, J. Chem. Phys. 148, 054902 (2018).

Model
-----
- The diffusing molecule cannot enter the lipid bilayer (short-range hard repulsion).
  Its reduced diffusion coefficient D/D0 in the network of water channels is
  interpolated from Brownian Dynamics simulations as a function of x = L/a
  (L: lipid length, a: lattice parameter); see ``tables.py``.
- A layer of water of thickness w next to the lipid heads is tightly bound to them
  and does not contribute to diffusion. Its fraction x_b of the total water rescales
  the result: D = (D/D0)(x) * (1 - x_b) * D0.
- The reference D0, measured in pure water at temperature T0, is rescaled to the
  experimental temperature Texp with the Stokes-Einstein relation (D ~ T / eta),
  using the Vogel equation for the viscosity of water.

Units: lengths in nm, water content in % weight, lipid density in g/cm^3,
temperatures in degrees Celsius, D0 and D in 10^-5 cm^2/s.
"""

import math
from dataclasses import dataclass

from .tables import TABLES

#: Minimal-surface parameters (A0: surface area of the unit cell in units of
#: a^2; chi: Euler-Poincare characteristic) of the bicontinuous cubic phases.
CUBIC_GEOMETRY = {
    "Pn3m": (1.919, -2),
    "Im3m": (2.345, -4),
    "Ia3d": (3.091, -8),
}

PHASES = ("Pn3m", "Im3m", "Ia3d")


@dataclass
class DiffusionResult:
    """Result of the calculation.

    bound_water_fraction:  x_b, # bound water molecules / total # water molecules
    D:                     effective diffusion coefficient (10^-5 cm^2/s)
    D0_at_Texp:            diffusion coefficient in pure water at Texp (10^-5 cm^2/s)
    reduced_D:             (D/D0) of the bare channel geometry, before bound-water correction
    lipid_length:          lipid length (nm)
    """

    phase: str
    bound_water_fraction: float
    D: float
    D0_at_Texp: float
    reduced_D: float
    lipid_length: float


def water_viscosity(temperature_c):
    """Viscosity of water (Pa s) from the Vogel equation; temperature in Celsius."""
    return 10 ** (-4.5318 - 220.57 / (149.39 - temperature_c - 273.15))


def rescale_D0(D0, T0, Texp):
    """Rescale D0 measured at T0 to Texp (Celsius) with the Stokes-Einstein relation."""
    return D0 * water_viscosity(T0) / (T0 + 273.15) * (Texp + 273.15) / water_viscosity(Texp)


def water_volume_fraction(water_weight, lipid_density):
    """Volume fraction of water from its weight percentage and the lipid density."""
    if not 0 < water_weight < 100:
        raise ValueError("Water content must be between 0 and 100 % weight")
    if lipid_density <= 0:
        raise ValueError("Lipid density must be positive")
    return 1 - 1 / (1 + lipid_density * water_weight / (100 - water_weight))


def _volume_fraction_within(phase, x):
    """Volume fraction within distance x*a from the minimal surface."""
    a0, chi = CUBIC_GEOMETRY[phase]
    return 2 * a0 * x + 4 / 3 * math.pi * chi * x ** 3


def cubic_lipid_length(lattice, phase, lipid_volume_fraction):
    """Lipid length (nm), solving phi_lipid = 2 A0 x + (4/3) pi chi x^3 for x = L/a."""
    a0, chi = CUBIC_GEOMETRY[phase]
    p = 3 * a0 / (2 * math.pi * chi)
    q = -3 * lipid_volume_fraction / (4 * math.pi * chi)
    arg = 3 * q / (2 * p) * math.sqrt(-3 / p)
    if not -1 <= arg <= 1:
        raise ValueError("No physical solution: lipid volume fraction %.3f too large for %s"
                         % (lipid_volume_fraction, phase))
    x = 2 * math.sqrt(-p / 3) * math.cos(math.acos(arg) / 3 - 2 * math.pi / 3)
    return x * lattice


def reduced_diffusion(phase, x):
    """D/D0 of a point-like tracer in the channels, linearly interpolated at x = L/a."""
    xs, ds = TABLES[phase]
    if x < xs[0] or x > xs[-1]:
        raise ValueError("L/a = %.4f is outside the simulated range [%g, %g] for %s"
                         % (x, xs[0], xs[-1], phase))
    if x == xs[0]:
        return ds[0]
    for i in range(1, len(xs)):
        if xs[i - 1] < x <= xs[i]:
            return ds[i - 1] + (ds[i] - ds[i - 1]) * (x - xs[i - 1]) / (xs[i] - xs[i - 1])


def diffusion(phase, lattice, water_weight, lipid_density, D0, T0, Texp, bound_layer):
    """Effective diffusion coefficient of a hydrophilic molecule in a cubic phase.

    Parameters
    ----------
    phase : str
        One of "Pn3m", "Im3m", "Ia3d".
    lattice : float
        Lattice parameter in nm, as measured by SAXS.
    water_weight : float
        Weight percentage of water used to prepare the cubic phase.
    lipid_density : float
        Density of the lipid in g/cm^3.
    D0 : float
        Diffusion coefficient of the molecule in pure water at T0 (10^-5 cm^2/s).
    T0, Texp : float
        Reference and experimental temperatures (Celsius).
    bound_layer : float
        Thickness w (nm) of the layer of bound water next to the lipid heads;
        typical values 0.6-1.2 nm (about two to four water molecules).
    """
    if phase not in PHASES:
        raise ValueError("Unknown phase %r; choose among %s" % (phase, ", ".join(PHASES)))
    if lattice <= 0:
        raise ValueError("Lattice parameter must be positive")
    if bound_layer < 0:
        raise ValueError("Bound-water thickness must be non-negative")

    phi_l = 1 - water_volume_fraction(water_weight, lipid_density)
    length = cubic_lipid_length(lattice, phase, phi_l)
    x = length / lattice
    reduced = reduced_diffusion(phase, x)

    a0, chi = CUBIC_GEOMETRY[phase]
    radius = math.sqrt(-a0 / (2 * math.pi * chi)) * lattice - length
    if bound_layer >= radius:
        raise ValueError("The bound-water layer (w = %g nm) fills the water channels "
                         "(channel radius %.3f nm)" % (bound_layer, radius))
    x_b = (_volume_fraction_within(phase, x + bound_layer / lattice) - phi_l) / (1 - phi_l)

    D0_exp = rescale_D0(D0, T0, Texp)
    return DiffusionResult(phase=phase, bound_water_fraction=x_b, D=reduced * (1 - x_b) * D0_exp,
                           D0_at_Texp=D0_exp, reduced_D=reduced, lipid_length=length)
