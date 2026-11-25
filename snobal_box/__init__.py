"""snobal-box: point-scale snobal for scenario testing and teaching.

An independent, pure-Python packaging of the snobal point energy and mass
balance snow model. The physics core is vendored from pysnobal 0.2.3 (CC0, a
work of the United States government); everything else is original. See
PROVENANCE.md.

Quickstart (config driven, the recommended way):

    import snobal_box as sb
    result = sb.run_config('examples/atwater_wy25_insitu.yaml')
    result.summary

or ``snobal-box run examples/atwater_wy25_insitu.yaml`` from a shell.
The older ``Site`` class still works for scripted use.
"""

__version__ = '0.2.0.dev0'

from snobal_box.io import (read_forcing, read_initial_snow,
                           read_measurement_heights, read_reference_output,
                           DATA_COLS, PRECIP_COLS, OUT_COLS)
from snobal_box.phase import (PHASE_TABLE, apply_phase_table,
                              shifted_phase_table, phase_table_frame,
                              repair_precip)
from snobal_box.perturb import perturb, snow_fraction
from snobal_box.model import (build_tstep_info, run, meltout, state_at,
                              swe_series)
from snobal_box.site import Site
from snobal_box.precip import (apply_precip_scheme, PHASE_SCHEMES,
                               DENSITY_SCHEMES, dewpoint_C)
from snobal_box.align import find_start
from snobal_box.metrics import season_metrics, compare, mass_balance
from snobal_box.config import load_config
from snobal_box.runner import (run_config, load_forcing, scenario_forcing,
                               run_forcing, site_state, RunResult)
from snobal_box.qc import qc_report
from snobal_box.diagnostics import (stats, find_offset, sensitivity_table,
                                    error_grid, plot_error_grid,
                                    COMPARE_PAIRS)

__all__ = [
    'run_config', 'load_config', 'load_forcing', 'scenario_forcing',
    'run_forcing', 'site_state', 'RunResult', 'apply_precip_scheme',
    'PHASE_SCHEMES', 'DENSITY_SCHEMES', 'dewpoint_C', 'find_start',
    'season_metrics', 'compare', 'mass_balance',
    '__version__', 'Site', 'run', 'perturb', 'qc_report', 'stats',
    'build_tstep_info', 'meltout', 'state_at', 'swe_series', 'snow_fraction',
    'find_offset', 'sensitivity_table', 'error_grid', 'plot_error_grid',
    'PHASE_TABLE', 'apply_phase_table', 'shifted_phase_table',
    'phase_table_frame', 'repair_precip', 'read_forcing', 'read_initial_snow',
    'read_measurement_heights', 'read_reference_output',
    'DATA_COLS', 'PRECIP_COLS', 'OUT_COLS', 'COMPARE_PAIRS',
]
