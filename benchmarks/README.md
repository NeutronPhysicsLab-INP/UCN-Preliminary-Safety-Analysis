# Benchmarks

Analytical validation cases for the AlSUN safety solver. Each notebook
runs against a closed-form reference and reports the deviation.

## Available benchmarks

| File | Reference | What it validates |
|---|---|---|
| [`blowdown_benchmark.ipynb`](blowdown_benchmark.ipynb) | Adiabatic, isentropic blowdown of an ideal gas through a choked orifice to vacuum. Closed-form: $m(t)/m_0 = (1 + \tfrac{\gamma-1}{2}\,\kappa\, t)^{-2/(\gamma-1)}$ | `alsun_safety.relief.flow.mdot_orifice` and the mass-balance step inside `step_burst_disk_absolute`. |

## Reproducing the validation

The notebooks ship with their output cells already populated, so the
results render directly on GitHub. To re-run:

```bash
pip install -e .                 # install the package + its deps
pip install jupyter              # if you don't already have Jupyter
jupyter notebook benchmarks/blowdown_benchmark.ipynb
# Run All Cells
```

Each notebook prints a `PASS:` or `REVIEW:` banner at the end based on
the maximum relative error against its analytical reference. The current
blowdown benchmark passes with $L^\infty$ relative error
$\sim 3 \times 10^{-5}$ for an explicit Euler integrator at
$\Delta t = 10$ ms.

## Adding new benchmarks

Each benchmark should:

1. State the closed-form reference solution at the top of the notebook,
   with full citation to the textbook or primary source.
2. Define the problem parameters in a single configuration cell.
3. Compute the analytical trajectory.
4. Compute the numerical trajectory using `alsun_safety` package
   functions (no reimplementation).
5. Plot the two on a shared time axis with a residual subplot.
6. Report $L^\infty$ and $L^2$ error norms.
7. Print a `PASS:` / `REVIEW:` banner against an explicit threshold.

Suggested next benchmarks:

- **Mariotte hoop stress** — verify
  `alsun_safety.stress.thin_wall.calculate_dynamic_stresses` against
  $\sigma_\theta = PR/t$ for a uniform-pressure thin-walled cylinder.
- **Bernoulli liquid discharge** — verify
  `alsun_safety.relief.flow.mdot_bernoulli` against
  $\dot m = C_d A \sqrt{2\rho \Delta P}$.
- **Detonation pressure-decay chain** — verify
  `alsun_safety.scenarios.detonation.run_detonation` against the
  closed-form sequential expansion $P_n / P_{n-1} = V_{n-1}/V_n$.

These have closed-form references and are short enough to add without
disrupting the notebook scaffolding.
