"""Command-line interface.

Single sample (defaults: glucose in a monolinolein Pn3m phase at 37 C, Antognini et al.):
    python -m mesophase_diffusion --phase Pn3m --lattice 8.57 --water 36 --density 1.05 \\
        --D0 0.67 --T0 25 --Texp 37 --w 0.8

Batch (CSV with columns phase,lattice,water,density,D0,T0,Texp,w; results written as CSV):
    python -m mesophase_diffusion --csv examples/samples.csv > results.csv
"""

import argparse
import csv
import sys

from .core import PHASES, diffusion

INPUTS = ["phase", "lattice", "water", "density", "D0", "T0", "Texp", "w"]
OUTPUTS = ["lipid_length_nm", "reduced_D", "D0_at_Texp", "bound_water_fraction", "D"]


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="mesophase-diffusion",
        description="Effective diffusion coefficient of a hydrophilic molecule in a lipidic "
                    "cubic phase (D, D0 in 10^-5 cm^2/s).")
    ap.add_argument("--phase", choices=PHASES, default="Pn3m", help="cubic phase (default Pn3m)")
    ap.add_argument("--lattice", type=float, default=8.57, help="lattice parameter, nm (8.57)")
    ap.add_argument("--water", type=float, default=36, help="water content, %% weight (36)")
    ap.add_argument("--density", type=float, default=1.05, help="lipid density, g/cm^3 (1.05)")
    ap.add_argument("--D0", type=float, default=0.67,
                    help="diffusion coefficient in pure water at T0, 10^-5 cm^2/s (0.67)")
    ap.add_argument("--T0", type=float, default=25, help="reference temperature, Celsius (25)")
    ap.add_argument("--Texp", type=float, default=37, help="experimental temperature, Celsius (37)")
    ap.add_argument("--w", type=float, default=0.8,
                    help="thickness of the bound-water layer, nm (0.8; typically 0.6-1.2)")
    ap.add_argument("--csv", metavar="FILE",
                    help="batch mode: CSV file with columns " + ",".join(INPUTS))
    args = ap.parse_args(argv)

    try:
        if args.csv:
            out = csv.writer(sys.stdout)
            out.writerow(INPUTS + OUTPUTS)
            with open(args.csv, newline="") as fh:
                for rec in csv.DictReader(fh):
                    vals = [float(rec[k]) for k in INPUTS[1:]]
                    res = diffusion(rec["phase"].strip(), *vals)
                    out.writerow([rec["phase"].strip()] + ["%.6g" % v for v in vals] +
                                 ["%.6g" % v for v in (res.lipid_length, res.reduced_D,
                                                       res.D0_at_Texp, res.bound_water_fraction,
                                                       res.D)])
        else:
            res = diffusion(args.phase, args.lattice, args.water, args.density,
                            args.D0, args.T0, args.Texp, args.w)
            print("Phase:                        %s" % res.phase)
            print("Lipid length (nm):            %.4f" % res.lipid_length)
            print("D/D0 of channel geometry:     %.4f" % res.reduced_D)
            print("D0 at Texp (1e-5 cm^2/s):     %.4f" % res.D0_at_Texp)
            print("Bound-water fraction x_b:     %.4f" % res.bound_water_fraction)
            print("D (1e-5 cm^2/s):              %.4f" % res.D)
    except ValueError as err:
        sys.exit("Error: %s" % err)


if __name__ == "__main__":
    main()
