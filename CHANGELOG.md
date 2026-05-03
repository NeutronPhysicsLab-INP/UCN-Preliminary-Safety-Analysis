# Changelog

All notable changes to the AlSUN safety solver will be documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and the project adheres to [Semantic Versioning](https://semver.org/).

## [1.0.0] — 2026-04-30

### Added
- Initial public release accompanying the manuscript
  *Preliminary Safety Assessment of Liquid Deuterium Premoderator
  Container for Ultra-Cold Neutron Source at WWR-K Reactor (AlSUN).*
- `alsun_safety` Python package:
  - `config` — single-source-of-truth `SimConfig` dataclass.
  - `geometry` — `R_curv` solver, free vacuum-volume calculation.
  - `eos` — NIST-tabulated equation-of-state handler with phase-detection.
  - `relief` — Bernoulli + API 520 orifice flow + rupture-disk dynamics.
  - `stress` — thin-walled vessel and two-wall LOVA stress analysis.
  - `scenarios` — cooling failure, LOVA (with and without helium gap),
    and bounding detonation case.
  - `uq` — multiprocessing-based Monte Carlo driver.
  - `plotting` — saturation curve, isotherm, MC-envelope helpers.
- Runnable entry-point scripts in `scripts/` for each accident scenario.
- Smoke tests in `tests/` exercising the no-NIST-data code paths.
