# Property data for the AlSUN safety solver

The solver requires two sets of NIST tabular property data for normal deuterium (D₂):

## Files expected in this directory

| File | Description |
|---|---|
| `fluid-2.cgi` | Saturation table (P, T, ρ, u, h, s, c_p, c_v, k for liquid and vapor along the saturation line) |
| `isotherms/T_<T>K.cgi` | Isothermal tables (one file per temperature, e.g., `T_30K.cgi`, `T_60K.cgi`, ..., `T_300K.cgi`) covering 20–300 K |

## How to obtain

Both are exported in plain tab-separated form from the NIST Chemistry WebBook:

<https://webbook.nist.gov/chemistry/fluid/>

Select fluid: **deuterium (R-704)**, choose temperature/pressure ranges, and "Tab-delimited" output. Save each one with the filename above.

The solver expects the columns NIST exports by default; the column-name parser in `alsun_safety/eos/property_db.py` is robust to the precise naming.

## Citation for the data

> Richardson, I. A., Leachman, J. W., Lemmon, E. W. (2014). *Fundamental Equation of State for Deuterium.*
> J. Phys. Chem. Ref. Data 43(1), 013103. DOI: 10.1063/1.4864752
