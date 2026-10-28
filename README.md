# MesophaseDiffusion

Estimates the effective diffusion coefficient of a hydrophilic molecule in the water channels of a lipidic bicontinuous cubic phase (Pn3m, Im3m, Ia3d) from the structural features of the phase. It follows the analysis of [Assenza and Mezzenga, *J. Chem. Phys.* 148, 054902 (2018)](https://doi.org/10.1063/1.5019224).

This is the Python version of the *Diffusion Coefficient in Water Channels* tool available online at
[www.salvatoreassenza.com/#tools](https://www.salvatoreassenza.com/#tools). It gives identical results, and it adds batch processing of whole data sets and use from your own scripts.

Pure Python (≥ 3.8), no dependencies.

## Installation

```bash
pip install git+https://github.com/saassenza/MesophaseDiffusion
```

or simply download the repository and run it from its folder (`python -m mesophase_diffusion ...`).

## Usage

### Command line

The default values correspond to glucose diffusing in a monolinolein Pn3m phase at 37 °C (Antognini et al., 2016):

```bash
mesophase-diffusion --phase Pn3m --lattice 8.57 --water 36 --density 1.05 \
    --D0 0.67 --T0 25 --Texp 37 --w 0.8
```

```
Phase:                        Pn3m
Lipid length (nm):            1.5052
D/D0 of channel geometry:     0.4721
D0 at Texp (1e-5 cm^2/s):     0.8993
Bound-water fraction x_b:     0.6480
D (1e-5 cm^2/s):              0.1494
```

Batch mode: a CSV file with columns `phase,lattice,water,density,D0,T0,Texp,w` (see `examples/samples.csv`):

```bash
mesophase-diffusion --csv examples/samples.csv > results.csv
```

If the package is not installed, replace `mesophase-diffusion` with `python -m mesophase_diffusion`.

### Python

```python
from mesophase_diffusion import diffusion

res = diffusion("Pn3m", lattice=8.57, water_weight=36, lipid_density=1.05,
                D0=0.67, T0=25, Texp=37, bound_layer=0.8)
print(res.D, res.bound_water_fraction)
```

## Inputs and units

| Input | Unit | Description |
|---|---|---|
| `phase` | – | `Pn3m`, `Im3m` or `Ia3d` |
| `lattice` | nm | lattice parameter *a*, from SAXS |
| `water` | % weight | water used to prepare the cubic phase |
| `density` | g/cm³ | density of the lipid (water is taken as 1 g/cm³) |
| `D0` | 10⁻⁵ cm²/s | diffusion coefficient of the molecule in pure water at `T0` |
| `T0` | °C | temperature at which `D0` is known |
| `Texp` | °C | temperature of the experiment |
| `w` | nm | thickness of the layer of bound water next to the lipid heads; good values are usually 0.6–1.2 nm (about two to four water molecules) |

Outputs: the fraction of bound water *x*<sub>b</sub> (bound water molecules / total water molecules) and the effective diffusion coefficient *D* (10⁻⁵ cm²/s).

## Model

1. **Structure.** The lipid length *L* follows from the lattice parameter and the lipid volume fraction φ<sub>L</sub>. It solves φ<sub>L</sub> = 2*A*<sub>0</sub>(*L*/*a*) + (4π/3)χ(*L*/*a*)³, with (*A*<sub>0</sub>, χ) = (1.919, −2), (2.345, −4) and (3.091, −8) for Pn3m, Im3m and Ia3d.
2. **Channel geometry.** The molecule cannot enter the lipid bilayer, so the network of water channels slows down its diffusion by a factor (*D*/*D*<sub>0</sub>)(*L*/*a*). This factor was obtained from Brownian Dynamics simulations (Assenza and Mezzenga, 2018) and is linearly interpolated from the tables in `mesophase_diffusion/tables.py`.
3. **Bound water.** A layer of thickness *w* next to the lipid heads is tightly bound to them and does not contribute to diffusion. Its fraction of the total water is

   $$x_b = \frac{2A_0\,\xi + \frac{4\pi}{3}\chi\,\xi^3 - \phi_L}{1-\phi_L}, \qquad \xi = \frac{L + w}{a}.$$

4. **Temperature.** *D*<sub>0</sub> is rescaled from *T*<sub>0</sub> to *T*<sub>exp</sub> with the Stokes–Einstein relation, *D* ∝ *T*/η(*T*). The viscosity of water comes from the Vogel equation, log<sub>10</sub> η = −4.5318 − 220.57 / (149.39 − *T*), with *T* in K and η in Pa s.
5. **Result:** *D* = (*D*/*D*<sub>0</sub>)(*L*/*a*) · (1 − *x*<sub>b</sub>) · *D*<sub>0</sub>(*T*<sub>exp</sub>).

### Range of validity

The simulation tables cover *L*/*a* ∈ [0, 0.3497] (Pn3m), [0, 0.2350] (Im3m) and [0, 0.2160] (Ia3d). The bound-water layer must also be thinner than the channel radius. Outside these limits the program stops with an explanatory error. The online tool instead returns `NaN` in the first case and meaningless numbers in the second; this is the only difference in behaviour.

## Tests

```bash
python -m unittest discover tests
```

The tests check that the results match the online tool.

## How to cite

If you use this tool in your research, please cite:

- S. Assenza and R. Mezzenga, *J. Chem. Phys.* **148**, 054902 (2018). [doi:10.1063/1.5019224](https://doi.org/10.1063/1.5019224)
- L. M. Antognini, S. Assenza, C. Speziale and R. Mezzenga, *J. Chem. Phys.* **145**, 084903 (2016). [doi:10.1063/1.4961224](https://doi.org/10.1063/1.4961224)

See also the companion tool [MesophaseStructure](https://github.com/saassenza/MesophaseStructure), which computes channel radii and lipid lengths for cubic, hexagonal and lamellar phases.

## License

MIT, see [LICENSE](LICENSE).

## Contact

Salvatore Assenza, Universidad Autónoma de Madrid (IFIMAC) – [www.salvatoreassenza.com](https://www.salvatoreassenza.com)
