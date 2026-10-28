"""Tests: reference values reproduce the online tool at https://www.salvatoreassenza.com/#tools.

Run with:  python -m unittest discover tests
"""

import subprocess
import sys
import unittest
from pathlib import Path

from mesophase_diffusion import diffusion, reduced_diffusion, rescale_D0, water_viscosity
from mesophase_diffusion.tables import TABLES

ROOT = Path(__file__).resolve().parents[1]

# (phase, lattice nm) -> (x_b, D in 1e-5 cm^2/s); 36 % water, density 1.05,
# D0 = 0.67 at 25 C, Texp = 37 C, w = 0.8 nm
REFERENCE = {
    ("Pn3m", 8.57): (0.6480346208587769, 0.1494122774361131),
    ("Im3m", 11.5): (0.5707510135460483, 0.14295831938636566),
    ("Ia3d", 13.0): (0.6954452403330647, 0.13566520592621675),
}


class TestDiffusion(unittest.TestCase):
    def test_reference_values(self):
        for (phase, lattice), (xb, D) in REFERENCE.items():
            res = diffusion(phase, lattice, 36, 1.05, 0.67, 25, 37, 0.8)
            self.assertAlmostEqual(res.bound_water_fraction, xb, places=12)
            self.assertAlmostEqual(res.D, D, places=12)

    def test_tables_consistent(self):
        for phase, (xs, ds) in TABLES.items():
            self.assertEqual(len(xs), len(ds), phase)
            self.assertTrue(all(b > a for a, b in zip(xs, xs[1:])), phase)

    def test_interpolation_hits_nodes(self):
        xs, ds = TABLES["Pn3m"]
        for i in (0, 10, len(xs) - 1):
            self.assertAlmostEqual(reduced_diffusion("Pn3m", xs[i]), ds[i], places=12)

    def test_viscosity_of_water(self):
        # ~0.89 mPa s at 25 C
        self.assertAlmostEqual(water_viscosity(25) * 1e3, 0.89, delta=0.01)
        self.assertAlmostEqual(rescale_D0(1.0, 30, 30), 1.0, places=12)

    def test_no_bound_layer(self):
        res = diffusion("Pn3m", 8.57, 36, 1.05, 0.67, 25, 25, 0.0)
        self.assertAlmostEqual(res.bound_water_fraction, 0.0, places=12)
        self.assertAlmostEqual(res.D, res.reduced_D * 0.67, places=12)

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):
            diffusion("Hex", 8.57, 36, 1.05, 0.67, 25, 37, 0.8)
        with self.assertRaises(ValueError):  # bound layer filling the channels
            diffusion("Pn3m", 8.57, 36, 1.05, 0.67, 25, 37, 5.0)
        with self.assertRaises(ValueError):  # L/a outside the simulated range
            diffusion("Ia3d", 5.0, 20, 1.05, 0.67, 25, 37, 0.8)

    def test_cli(self):
        out = subprocess.run([sys.executable, "-m", "mesophase_diffusion"], cwd=ROOT,
                             capture_output=True, text=True, check=True).stdout
        self.assertIn("0.1494", out)

    def test_cli_csv(self):
        out = subprocess.run([sys.executable, "-m", "mesophase_diffusion", "--csv",
                              "examples/samples.csv"],
                             cwd=ROOT, capture_output=True, text=True, check=True).stdout
        self.assertEqual(len(out.strip().splitlines()), 1 + 3)


if __name__ == "__main__":
    unittest.main()
