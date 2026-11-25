# snobal-box

[![tests](https://github.com/MacSkiles/snobal-in-a-box/actions/workflows/tests.yml/badge.svg)](https://github.com/MacSkiles/snobal-in-a-box/actions/workflows/tests.yml) [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/MacSkiles/snobal-in-a-box/blob/main/notebooks/snobal_box_quickstart.ipynb)

snobal-box runs the snobal point snow energy and mass balance model
(Marks et al., 1999) from a single YAML configuration file, using either the
bundled example forcings or user-supplied hourly data. It is pure Python and
requires no compiled code or gridded setup.

```
pip install git+https://github.com/MacSkiles/snobal-in-a-box
snobal-box init my_site                 # template configuration and forcing CSV
snobal-box check my_site/my_site.yaml   # read and quality-check the forcing
snobal-box run   my_site/my_site.yaml   # run baseline and scenarios
```

From Python:

```python
import snobal_box as sb
res = sb.run_config('examples/atwater_wy25_insitu.yaml')
res.summary                 # one row per scenario
res.outputs['baseline']     # hourly model output
```

The notebook `notebooks/snobal_box_quickstart.ipynb` runs in Google Colab
and covers both the bundled example and user-uploaded forcing.

## Inputs

**Forcing** is hourly, in one of two formats:

| Format | Description | Start time |
|---|---|---|
| `csv` | A single file with a time column. Column names are mapped by a preset (`snobal`, `hrrr_extract`) or a `columns:` mapping. | Taken from the timestamps |
| `input` | The snobal `.input` file pair (data and precipitation files). | Set in the configuration: `start` or `water_year` |

Required variables are net solar radiation `S_n`, incoming longwave
radiation `I_lw`, air temperature `T_a`, vapor pressure `e_a`, wind speed `u`,
soil temperature `T_g` and precipitation mass `m_pp`. The variables
`percent_snow`, `rho_snow` and `T_pp` are optional; when absent they are
derived by the selected precipitation scheme.

**Site** parameters are elevation, measurement heights (`z_u`, `z_t`, `z_0`,
`z_g`) and initial snow conditions. They are given in the configuration, or
read from snobal `inheight` and `snow.properties` `.input` files.

## Start time

snobal `.input` files contain an hourly counter rather than timestamps, so the
start time must be specified:

- `water_year: 2025` sets the first record to 2024-10-01 00:00 in `timezone`.
- `start: '2024-10-01 00:00'` sets it explicitly.

Each run writes `config_resolved.yaml`, recording the start time used.

All times are timezone-aware. `timezone` is the zone of the input times and
`output_timezone` the zone of the output (for example, HRRR forcing in UTC
reported in MST as `Etc/GMT+7`).

## Precipitation phase and new snow density

snobal requires the snow fraction and new snow density of precipitation as
inputs. The `precip:` block selects how they are determined:

```yaml
precip:
  scale: 1.0
  phase: table              # auto | file | table | threshold | linear | dai2008
  phase_params: {}          # threshold: {T_threshold: 0}; linear: {T_snow: -1, T_rain: 3}
  phase_temperature: T_pp   # T_pp | T_a | dewpoint
  density: table            # auto | file | table | constant | hedstrom_pomeroy
  density_params: {}
  threshold_shift: 0.0      # degrees C; positive values increase the snow fraction
```

`file` uses the values in the forcing. `auto` (the default) uses `file` when
those columns are present and `table` otherwise. `table` is a temperature
lookup table consistent with the Atwater precipitation files. Additional
schemes can be registered in `snobal_box/precip.py`.

## Scenarios

Each scenario applies changes to the baseline forcing:

```yaml
scenarios:
  +2C:               {T_offset: 2.0}                         # warming, with phase recomputed
  +2C energy only:   {T_offset: 2.0, rebuild_phase: false}   # warming, baseline phase retained
  +2C phase only:    {threshold_shift: -2.0}                 # phase as under +2 C, no warming
  dai2008 phase:     {precip: {phase: dai2008, density: hedstrom_pomeroy}}
  precip +50%:       {precip_scale: 1.5}
```

Available parameters: `T_offset`, `solar_scale`, `lw_offset`, `wind_scale`,
`vapor_scale`, `precip_scale`, `rho_snow_scale`, `threshold_shift`,
`rebuild_phase` and `precip`.

## Outputs

Each run writes to `output/<name>/`:

- `summary.csv`, one row per scenario: peak SWE and date, half-melt date,
  melt-out date, depletion duration and rate, peak depth, snowfall, rain on
  snow, rain on bare ground, snow water input, sublimation, mass balance
  closure residual, snow fraction, number of neutral-stability substitutions,
  and peak depth relative to wind measurement height. When observations are
  supplied, bias, RMSE, correlation and timing errors for depth and SWE.
- `output_<scenario>.csv`, hourly model state and fluxes.
- `swe_depth.png`, `water.png` and `energy.png`.
- `config_resolved.yaml`, the full configuration used for the run.

`snobal-box run --combined all.csv` combines the summaries of several
configurations.

## Examples

The configurations in `examples/` run against the data in `data/`:

| Configuration | Forcing |
|---|---|
| `atwater_wy25_insitu.yaml` | Atwater Study Plot station, WY2025, `.input` files |
| `atwater_wy24_insitu.yaml` | Atwater Study Plot station, WY2024, `.input` files |
| `atwater_wy24_hrrr.yaml` | HRRR extracted at the plot, WY2024, CSV, with four phase and density schemes |

WY2024 has both station and HRRR forcing, so the two WY2024 examples can be
compared directly.

## Provenance and license

Code is released under the MIT License (`LICENSE`). The physics core is
PySnobal 0.2.3 (CC0, a work of the US government); see `PROVENANCE.md`.
Data in `data/` are released under CC-BY-4.0 (`data/LICENSE`). A Zenodo DOI
will be added when the forcings are archived.

Sites: Atwater Study Plot (UT). Senator Beck and Swamp Angel Study Plots (CO)
will follow.

## Reference

Marks, D., Domingo, J., Susong, D., Link, T., and Garen, D. (1999). A
spatially distributed energy balance snowmelt model for application in
mountain basins. Hydrological Processes, 13, 1935–1959.
