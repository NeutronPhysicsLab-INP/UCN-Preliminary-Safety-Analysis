# AlSUN Safety Solver

Python solver accompanying the manuscript

> *Preliminary Safety Assessment of Liquid Deuterium Premoderator Container for Ultra-Cold Neutron Source at WWR-K Reactor (AlSUN)*

This repository contains the analytical / numerical solver used in the manuscript to evaluate the structural and thermal-hydraulic response of the liquid deuterium (LD₂) premoderator container under three categories of postulated accidents:

1. **Cooling failure** — loss of helium refrigeration with passive blowdown into a ballast tank.
2. **Loss-of-Vacuum Accident (LOVA)** — air ingress through the vacuum jacket with cryopumping deposition on cold surfaces; both with and without the helium coolant gap.
3. **Detonation** — beyond-design-basis residual-risk bound, sequential pressure decay through the layered shell geometry.

A Monte Carlo uncertainty quantification (UQ) driver wraps each scenario.

## Redacted release

A small set of input parameters derived from the WWR-K reactor's classified neutronic and structural data have been replaced by non-sensitive placeholder values in `alsun_safety/config.py`. These are flagged `# REDACTED` in the source. The full physics framework — equation-of-state handler, lumped-parameter transient solver, rupture-disk and ballast-tank dynamics, LOVA cryopumping model, structural stress chain, and Monte Carlo driver — is provided in source.

## Repository layout

```
github/
├── README.md                       this file
├── requirements.txt                Python dependencies
├── alsun_safety/                   the package
│   ├── config.py                   physical constants, geometry, scenario parameters
│   ├── eos/                        NIST-tabulated equation-of-state handler
│   │   ├── property_db.py          PropertyDB class, NIST loader
│   │   └── solver.py               state_from_T_rho, invert_T_from_u
│   ├── relief/                     rupture-disk + orifice flow
│   │   ├── flow.py                 Bernoulli + API 520 compressible flow
│   │   └── rupture_disk.py         step_burst_disk_absolute
│   ├── stress/                     thin-walled vessel stress
│   │   └── thin_wall.py            calculate_dynamic_stresses, calculate_2wall_stresses
│   ├── scenarios/                  per-scenario simulators
│   │   ├── cooling_failure.py      run_simulation_absolute_UQ
│   │   ├── lova_helium.py          run_simulation_Helium_Gap (with He gap)
│   │   ├── lova_no_helium.py       run_simulation_LOVA_REAL (single shell)
│   │   └── detonation.py           sequential layer-by-layer stress chain
│   ├── uq/
│   │   └── monte_carlo.py          multiprocessing MC driver
│   └── plotting/                   reusable plot helpers
├── scripts/                        runnable entry points
│   ├── run_cooling_failure.py
│   ├── run_lova_with_helium.py
│   ├── run_lova_without_helium.py
│   ├── run_detonation.py
│   └── air_ingress.py              R4-11 vacuum free volume + ingress rate
└── data/                           NIST tabular data (place .cgi files here)
    └── README.md                   instructions for obtaining the NIST data
```

## Installation

```bash
git clone <repo-url>
cd alsun-safety-solver
pip install -r requirements.txt
```

NIST saturation and isothermal property tables for deuterium are required (free download from <https://webbook.nist.gov/chemistry/fluid/>). Place them in `data/` — see `data/README.md` for the expected filenames.

## Quick start

```python
from alsun_safety.config import default_config
from alsun_safety.scenarios.cooling_failure import run_cooling_failure
from alsun_safety.eos.property_db import PropertyDB

cfg = default_config()
db  = PropertyDB.from_files(cfg.paths.saturation_file, cfg.paths.isotherm_dir)
result = run_cooling_failure(db, cfg)
```

Or run the prepared scripts directly:

```bash
python scripts/run_cooling_failure.py
python scripts/run_lova_with_helium.py
python scripts/run_detonation.py
python scripts/air_ingress.py
```

## Notes on reproducibility

- Numerical results in the manuscript figures were generated with the parameters in `config.py` set to the values listed in Table 1 of the manuscript.
- The Monte Carlo driver uses `multiprocessing.Pool` and reseeds each worker. Ensemble sizes default to N=100 per the manuscript; raise to ~1000 for tighter percentile estimates.
- The detonation analysis is closed-form (no time integration) and runs in seconds; the LOVA and cooling-failure simulations take ~30 s per realization on a desktop CPU.

## License

Add your preferred license. (MIT recommended for code; CC-BY for data/figures.)

## Citation

If you use this code, please cite the manuscript:

> Almukhametov, A.; Turlybekuly, K.; Shaimerdenov, A.; et al.
> *Preliminary Safety Assessment of Liquid Deuterium Premoderator Container for Ultra-Cold Neutron Source at WWR-K Reactor (AlSUN).*
> [Journal name], vol., pages, year. DOI:
