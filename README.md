# snobal-box

[![tests](https://github.com/MacSkiles/snobal-in-a-box/actions/workflows/tests.yml/badge.svg)](https://github.com/MacSkiles/snobal-in-a-box/actions/workflows/tests.yml) [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/MacSkiles/snobal-in-a-box/blob/main/notebooks/snobal_box_quickstart.ipynb)

Run the **snobal** point snow energy and mass balance model with your own
forcing or the bundled examples, from one YAML file. No compiled code, no
gridded setup: `pip install`, write a config, run.

```
pip install git+https://github.com/MacSkiles/snobal-in-a-box
snobal-box init my_site          # template config + template forcing CSV
snobal-box check my_site/my_site.yaml   # read and QC the forcing, no run
snobal-box run   my_site/my_site.yaml   # baseline + scenarios, CSV + figures
```

or in Python / a notebook:

```python
import snobal_box as sb
res = sb.run_config('examples/atwater_wy22_insitu.yaml')
res.summary                 # one row per scenario
res.outputs['baseline']     # hourly model output
```

A Colab-ready notebook is in `notebooks/snobal_box_quickstart.ipynb`: run the
example, or upload your own CSV and edit one cell.

## What goes in

**Forcing**, hourly, in one of two formats:

| format | what | start time |
|---|---|---|
| `csv` | one file with a time column. The upload format. Column names are mapped by a preset (`snobal`, `hrrr_extract`) or a `columns:` dict | from the timestamps |
| `ipw` | the classic snobal `.input` pair (data + precip files) | **you supply it**: `start`, `water_year`, or `start: auto` + `align_with` |

Required variables: net solar `S_n`, incoming longwave `I_lw`, air
temperature `T_a`, vapor pressure `e_a`, wind `u`, soil temperature `T_g`,
precipitation `m_pp`. Optional: `percent_snow`, `rho_snow`, `T_pp`; if
absent, the precipitation scheme makes them.

**Site**: elevation, instrument heights (`z_u`, `z_t`, `z_0`, `z_g`) and the
snow on the ground at the start. All in the config; IPW `inheight` and
`snow.properties` files are also accepted.

## Start times, made explicit

`.input` files carry an hourly counter, not times, so the start is an
assertion that nothing used to check. Now:

- `water_year: 2022` means 2021-10-01 00:00 in `timezone`.
- `align_with: <timestamped CSV>` slides the record's air temperature along
  the reference's (diurnal anomalies, hour resolution) and reports the best
  start and its correlation. With an explicit start it warns at 2 h or more
  of disagreement; with `start: auto` it uses the aligned start.
- Every run writes `config_resolved.yaml` with the start actually used and the
  alignment result.

Times are timezone-aware throughout. `timezone` is the zone of the input
times; `output_timezone` is the zone of the output (e.g. HRRR in UTC,
reported in MST as `Etc/GMT+7`).

## Precipitation phase and new snow density, made swappable

snobal does not decide rain versus snow; it is told. The `precip:` block
chooses who tells it:

```yaml
precip:
  scale: 1.0
  phase: table              # auto | file | table | threshold | linear | dai2008
  phase_params: {}          # threshold: {T_threshold: 0}; linear: {T_snow: -1, T_rain: 3}
  phase_temperature: T_pp   # T_pp | T_a | dewpoint
  density: table            # auto | file | table | constant | hedstrom_pomeroy
  density_params: {}
  threshold_shift: 0.0      # C, positive = boundary warmer = MORE snow
```

`file` uses whatever the forcing already contains; `auto` (default) is `file`
when the columns exist and `table` when they do not. `table` is the lookup
table reverse-engineered from the Atwater precip files (it reproduces their
snow fraction: 97% WY22, 94% WY24). Adding a scheme is one function and one
dictionary entry in `snobal_box/precip.py`.

## Scenarios

Each scenario is the baseline plus changes:

```yaml
scenarios:
  +2C:               {T_offset: 2.0}                         # warmer air AND more rain
  +2C energy only:   {T_offset: 2.0, rebuild_phase: false}   # warmer air, same rain/snow split
  +2C phase only:    {threshold_shift: -2.0}                 # same air, split as if 2 C warmer
  dai2008 phase:     {precip: {phase: dai2008, density: hedstrom_pomeroy}}
  HRRR +50%:         {precip_scale: 1.5}
```

Knobs: `T_offset`, `solar_scale`, `lw_offset`, `wind_scale`, `vapor_scale`,
`precip_scale`, `rho_snow_scale`, `threshold_shift`, `rebuild_phase`,
`precip: {...}`.

## What comes out

Per run, in `output/<name>/`:

- `summary.csv`, one row per scenario: peak SWE and date, half-melt date,
  melt-out date, depletion days and rate, peak depth, snowfall, rain on snow,
  rain on bare ground, snow water input, sublimation, **mass closure
  residual** (should be a few mm at most), snow fraction, neutral-stability
  substitutions, depth as a fraction of the wind sensor height, and, with
  observations, bias/RMSE/r and timing errors for depth and SWE.
- `output_<scenario>.csv`, hourly model state and fluxes.
- `swe_depth.png`, `water.png`, `energy.png`.
- `config_resolved.yaml`, everything needed to reproduce the run.

`--combined all.csv` stacks the summaries of several configs.

## Examples

`examples/` runs out of the box against `data/`:

| config | forcing |
|---|---|
| `atwater_wy22_insitu.yaml` | Atwater Study Plot station, WY2022, IPW |
| `atwater_wy24_insitu.yaml` | same, WY2024, using the corrected precip file (`_v2`: the original's precipitation is 24 h late) |
| `atwater_wy22_hrrr.yaml` | HRRR extracted at the plot, WY2022, CSV, with four phase/density schemes |

## Things that will bite you

- **A mid-winter start from bare ground is silently wrong.** Set
  `initial_snow.depth_m` whenever the record does not start before the first
  snow.
- **Depth cannot separate precipitation from new snow density.** Tuning
  `precip_scale` against depth assumes the density scheme is right.
- **The wind sensor is a ceiling.** snobal measures wind height above the
  snow surface; the summary reports peak depth over `z_u`.
- **No albedo model.** `S_n` is net solar, so it already contains the albedo
  of the pack that existed. Changing snowfall does not change brightness.
- **Neutral-stability substitutions.** When the turbulent-transfer iteration
  fails to converge, neutral fluxes are used for that step and counted.

## Provenance and license

Code: MIT (`LICENSE`). The physics core is PySnobal 0.2.3 (CC0, US
government work); see `PROVENANCE.md`.
Data in `data/`: CC-BY-4.0 (`data/LICENSE`). A Zenodo DOI will be added when
the forcings are archived.

Sites: Atwater Study Plot (UT) now; Senator Beck and Swamp Angel Study Plots
(CO) to follow.
