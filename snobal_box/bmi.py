"""
BMI 2.0 wrapper for point snobal.

Implements the CSDMS Basic Model Interface against the same physics core that
driver.py uses. The wrapper owns no physics: it is a translation layer between
BMI's set/update/get protocol and snobal's two-endpoint timestep call.

    from bmi_snobal import BmiSnobal
    m = BmiSnobal()
    m.initialize('atwater.yaml')
    while m.get_current_time() < m.get_end_time():
        m.set_value('land_surface_air__temperature', np.array([270.0]))
        ...
        m.update()
        swe = m.get_value('snowpack__liquid_equivalent_depth',
                          np.empty(1))
    m.finalize()

NOTE ON VARIABLE NAMES
The names below follow CSDMS Standard Name *style* and reuse the NextGen
forcing names where they exist, but they have NOT been validated against the
CSDMS Standard Names registry. Check them there before publishing, and expect
to change some. The mapping lives in one dict so that is a cheap edit.

NOTE ON TIMESTEP SEMANTICS
snobal advances across an interval using forcing at BOTH endpoints, and
interpolates between them when it subdivides. BMI has no way to express that,
so this wrapper buffers: values set via set_value() describe the END of the
interval that the next update() will cross. initialize() primes the start
endpoint from the config. This is documented here because a coupler that
assumes set_value() means "forcing during the coming step" will be off by one.
"""

from __future__ import annotations

import os
import warnings

import numpy as np
import pandas as pd
import yaml
from bmipy import Bmi

from snobal_box.site import Site
from snobal_box.model import build_tstep_info
from snobal_box._pysnobal.point import Snobal
from snobal_box._pysnobal.point.input import InputData

SECONDS_PER_HOUR = 3600.0


# --------------------------------------------------------------- variable map
# bmi name -> (internal key, units, role)
# `role` is 'forcing' for snobal input records, 'state' for SnowState attrs.
VARIABLES = {
    # inputs
    'land_surface_radiation~net~shortwave__energy_flux':
        ('S_n', 'W m-2', 'forcing'),
    'land_surface_radiation~incoming~longwave__energy_flux':
        ('I_lw', 'W m-2', 'forcing'),
    'land_surface_air__temperature':
        ('T_a', 'K', 'forcing'),
    'land_surface_air_water~vapor__partial_pressure':
        ('e_a', 'Pa', 'forcing'),
    'land_surface_wind__speed':
        ('u', 'm s-1', 'forcing'),
    'soil_surface__temperature':
        ('T_g', 'K', 'forcing'),
    'atmosphere_water__precipitation_leq_mass':
        ('m_pp', 'kg m-2', 'forcing'),
    'atmosphere_water__snowfall_mass_fraction':
        ('percent_snow', '1', 'forcing'),
    'snow~new__density':
        ('rho_snow', 'kg m-3', 'forcing'),
    'atmosphere_water__precipitation_temperature':
        ('T_pp', 'K', 'forcing'),

    # outputs
    'snowpack__liquid_equivalent_depth': ('m_s', 'kg m-2', 'state'),
    'snowpack__depth': ('z_s', 'm', 'state'),
    'snowpack__bulk_density': ('rho', 'kg m-3', 'state'),
    'snowpack__liquid_water_mass': ('h2o', 'kg m-2', 'state'),
    'snowpack__temperature': ('T_s', 'K', 'state'),
    'snowpack_surface__temperature': ('T_s_0', 'K', 'state'),
    'snowpack__cold_content': ('cc_s', 'J m-2', 'state'),
    'snowpack__melt_mass': ('melt', 'kg m-2', 'state'),
    'snowpack__runoff_mass': ('ro_predict', 'kg m-2', 'state'),
    'snowpack__sublimation_mass': ('E_s', 'kg m-2', 'state'),
    'snowpack__net_energy_flux': ('delta_Q_bar', 'W m-2', 'state'),
    'land_surface_net_radiation__energy_flux': ('R_n_bar', 'W m-2', 'state'),
    'land_surface_sensible_heat__energy_flux': ('H_bar', 'W m-2', 'state'),
    'land_surface_latent_heat__energy_flux': ('L_v_E_bar', 'W m-2', 'state'),
}

INPUT_VARS = tuple(n for n, v in VARIABLES.items() if v[2] == 'forcing')
OUTPUT_VARS = tuple(n for n, v in VARIABLES.items() if v[2] == 'state')


class BmiSnobal(Bmi):
    """BMI 2.0 interface to the point snobal energy and mass balance model."""

    _name = 'snobal (point energy and mass balance snow model)'
    _bmi_version = '2.0'

    def __init__(self):
        self._model = None
        self._buffers = {}      # bmi name -> 1-element float64 array
        self._next = {}         # forcing for the end of the coming interval
        self._current = {}      # forcing at the current model time
        self._time = 0.0        # seconds since start
        self._end_time = 0.0
        self._dt = SECONDS_PER_HOUR
        self._start_datetime = None
        self._forcing = None    # optional: drive from a file
        self._index = 0

    # ------------------------------------------------------------- bookkeeping

    def _check(self, name):
        if name not in VARIABLES:
            raise ValueError(f'unknown variable {name!r}')
        return VARIABLES[name]

    def _sync_outputs(self):
        state = self._model.snow_state
        for name in OUTPUT_VARS:
            key = VARIABLES[name][0]
            self._buffers[name][0] = float(getattr(state, key, np.nan))

    # ----------------------------------------------------------------- control

    def initialize(self, config_file=None):
        """Read the config, build the model, prime the first forcing endpoint.

        The config is YAML. A minimal one:

            site:
              root: data
              data_file: Atwater_snobal_Data.input
              precip_file: Atwater_snobal_PrecipFile.input
              snow_file: Atwater_snow_properties.input
              height_file: Atwater_inheight_properties_1.input
              elevation: 2600
              start_date: '2023-10-01'
              start_record: 686
            drive_from_file: true

        With `drive_from_file: true` the wrapper reads forcing from the files
        itself, so update() needs no set_value() calls. That is the mode to use
        for standalone runs and regression tests. Set it false when a coupler
        supplies forcing.
        """
        if config_file is None or not os.path.exists(config_file):
            raise FileNotFoundError(f'config file not found: {config_file!r}')
        with open(config_file) as fh:
            cfg = yaml.safe_load(fh)

        site = Site(name=cfg['site'].get('name', 'snobal'),
                           **{k: v for k, v in cfg['site'].items()
                              if k != 'name'})
        forcing, snow_prop, meas_heights, _ = site.load(verbose=False)

        self._drive_from_file = bool(cfg.get('drive_from_file', True))
        self._forcing = forcing
        self._start_datetime = forcing.index[0]
        self._dt = (forcing.index[1] - forcing.index[0]).total_seconds()
        self._end_time = (len(forcing) - 1) * self._dt
        self._time = 0.0
        self._index = 0

        tstep_info = build_tstep_info(data_tstep=int(self._dt / 60))
        params = {
            'start_date': self._start_datetime,
            'elevation': site.elevation,
            'data_tstep': self._dt / 60,
            'max_h2o_vol': site.max_h2o_vol,
            'max_z_s_0': site.max_z_s_0,
            'relative_heights': False,
        }
        self._model = Snobal(params, tstep_info, snow_prop, meas_heights)

        # one buffer per variable, so get_value_ptr can hand out a reference
        self._buffers = {name: np.zeros(1, dtype=np.float64)
                         for name in VARIABLES}

        # prime the start endpoint from the first forcing record
        first = forcing.iloc[0].to_dict()
        for name in INPUT_VARS:
            self._buffers[name][0] = first[VARIABLES[name][0]]
        self._current = dict(first)
        self._sync_outputs()

    def update(self):
        """Advance one data timestep.

        Forcing currently set describes the END of the interval being crossed.
        In drive_from_file mode the wrapper supplies it instead.
        """
        if self._model is None:
            raise RuntimeError('initialize() first')
        if self._index + 1 >= len(self._forcing):
            raise RuntimeError('model is already at its end time')

        if self._drive_from_file:
            end = self._forcing.iloc[self._index + 1].to_dict()
        else:
            end = {VARIABLES[n][0]: float(self._buffers[n][0])
                   for n in INPUT_VARS}

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            self._model.do_data_tstep(InputData(self._current),
                                      InputData(end))

        self._current = dict(end)
        self._index += 1
        self._time += self._dt
        self._sync_outputs()

    def update_until(self, time):
        if time < self._time:
            raise ValueError('cannot run backwards')
        while self._time + self._dt <= time + 1e-9:
            self.update()

    def finalize(self):
        self._model = None
        self._forcing = None
        self._buffers = {}

    # -------------------------------------------------------------------- info

    def get_component_name(self):
        return self._name

    def get_bmi_version(self):
        return self._bmi_version

    def get_input_item_count(self):
        return len(INPUT_VARS)

    def get_output_item_count(self):
        return len(OUTPUT_VARS)

    def get_input_var_names(self):
        return INPUT_VARS

    def get_output_var_names(self):
        return OUTPUT_VARS

    # --------------------------------------------------------------- variables

    def get_var_grid(self, name):
        self._check(name)
        return 0

    def get_var_type(self, name):
        self._check(name)
        return 'float64'

    def get_var_units(self, name):
        return self._check(name)[1]

    def get_var_itemsize(self, name):
        self._check(name)
        return np.dtype('float64').itemsize

    def get_var_nbytes(self, name):
        return self.get_var_itemsize(name) * self.get_grid_size(0)

    def get_var_location(self, name):
        self._check(name)
        return 'node'

    # -------------------------------------------------------------------- time

    def get_current_time(self):
        return self._time

    def get_start_time(self):
        return 0.0

    def get_end_time(self):
        return self._end_time

    def get_time_units(self):
        return 's'

    def get_time_step(self):
        return self._dt

    def get_current_datetime(self):
        """Not BMI. Convenience for callers that want a real timestamp."""
        return self._start_datetime + pd.to_timedelta(self._time, unit='s')

    # --------------------------------------------------------- getters/setters

    def get_value(self, name, dest):
        dest[:] = self._buffers[self._check(name) and name]
        return dest

    def get_value_ptr(self, name):
        self._check(name)
        return self._buffers[name]

    def get_value_at_indices(self, name, dest, inds):
        dest[:] = self.get_value_ptr(name)[inds]
        return dest

    def set_value(self, name, values):
        self._check(name)
        if name in OUTPUT_VARS:
            raise ValueError(f'{name!r} is an output variable')
        self._buffers[name][:] = np.asarray(values, dtype=np.float64).ravel()

    def set_value_at_indices(self, name, inds, src):
        self._check(name)
        if name in OUTPUT_VARS:
            raise ValueError(f'{name!r} is an output variable')
        self._buffers[name][inds] = src

    # -------------------------------------------------------------------- grid
    # A point model. One grid, id 0, type 'scalar', rank 0, size 1.

    def get_grid_type(self, grid):
        self._grid_check(grid)
        return 'scalar'

    def get_grid_rank(self, grid):
        self._grid_check(grid)
        return 0

    def get_grid_size(self, grid):
        self._grid_check(grid)
        return 1

    def get_grid_node_count(self, grid):
        self._grid_check(grid)
        return 1

    @staticmethod
    def _grid_check(grid):
        if grid != 0:
            raise ValueError(f'no grid {grid}')

    # Not meaningful for a scalar grid. The spec allows raising here.
    def get_grid_shape(self, grid, shape):
        raise NotImplementedError('scalar grid has no shape')

    def get_grid_spacing(self, grid, spacing):
        raise NotImplementedError('scalar grid has no spacing')

    def get_grid_origin(self, grid, origin):
        raise NotImplementedError('scalar grid has no origin')

    def get_grid_x(self, grid, x):
        raise NotImplementedError('scalar grid has no coordinates')

    def get_grid_y(self, grid, y):
        raise NotImplementedError('scalar grid has no coordinates')

    def get_grid_z(self, grid, z):
        raise NotImplementedError('scalar grid has no coordinates')

    def get_grid_edge_count(self, grid):
        raise NotImplementedError('scalar grid is not unstructured')

    def get_grid_face_count(self, grid):
        raise NotImplementedError('scalar grid is not unstructured')

    def get_grid_edge_nodes(self, grid, edge_nodes):
        raise NotImplementedError('scalar grid is not unstructured')

    def get_grid_face_edges(self, grid, face_edges):
        raise NotImplementedError('scalar grid is not unstructured')

    def get_grid_face_nodes(self, grid, face_nodes):
        raise NotImplementedError('scalar grid is not unstructured')

    def get_grid_nodes_per_face(self, grid, nodes_per_face):
        raise NotImplementedError('scalar grid is not unstructured')
