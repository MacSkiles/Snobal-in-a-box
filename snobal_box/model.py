"""Timestep configuration and the model run loop."""

import numpy as np
import pandas as pd

from snobal_box._pysnobal.point import Snobal
from snobal_box._pysnobal.point.input import InputData
from snobal_box._pysnobal.core.constants import (
    DATA_TSTEP, NORMAL_TSTEP, MEDIUM_TSTEP, SMALL_TSTEP, FREEZE)

__all__ = ['build_tstep_info', 'run', 'meltout', 'state_at', 'swe_series']

def build_tstep_info(data_tstep=60, norm=60, medium=15, small=1,
                     norm_threshold=60.0, medium_threshold=10.0,
                     small_threshold=1.0):
    """Replicates PySnobal.parse_time_step() without inicheck."""
    tstep_info = [{
        'level': i, 'output': False, 'threshold': None, 'time_step': None,
        'time_step_timedelta': None, 'intervals': None,
        'intervals_per_timestep': 1,
    } for i in range(4)]

    tstep_info[DATA_TSTEP].update({
        'time_step': data_tstep * 60,
        'time_step_timedelta': pd.to_timedelta(data_tstep, unit='min'),
    })
    tstep_info[NORMAL_TSTEP].update({
        'time_step': norm * 60,
        'time_step_timedelta': pd.to_timedelta(norm, unit='min'),
        'intervals': int(data_tstep / norm),
        'intervals_per_timestep': int(data_tstep / norm),
        'threshold': norm_threshold,
        'output': True,
    })
    tstep_info[MEDIUM_TSTEP].update({
        'time_step': medium * 60,
        'time_step_timedelta': pd.to_timedelta(medium, unit='min'),
        'intervals': int(norm / medium),
        'intervals_per_timestep': int(norm / medium),
        'threshold': medium_threshold,
    })
    tstep_info[SMALL_TSTEP].update({
        'time_step': small * 60,
        'time_step_timedelta': pd.to_timedelta(small, unit='min'),
        'intervals': int(medium / small),
        'intervals_per_timestep': int(norm / small),
        'threshold': small_threshold,
    })
    return tstep_info


def run(forcing, snow_prop, meas_heights, elevation,
        max_h2o_vol=0.01, max_z_s_0=0.25, tstep_info=None, verbose=False):
    """Run the point model over a forcing DataFrame. Returns a DataFrame."""
    if tstep_info is None:
        tstep_info = build_tstep_info()

    params = {
        'start_date': forcing.index[0],
        'elevation': elevation,
        'data_tstep': tstep_info[DATA_TSTEP]['time_step'] / 60,
        'max_h2o_vol': max_h2o_vol,
        'max_z_s_0': max_z_s_0,
        'relative_heights': False,
    }

    model = Snobal(params, tstep_info, snow_prop, meas_heights)

    rows = forcing.iterrows()
    _, input1 = next(rows)
    for index, input2 in rows:
        model.do_data_tstep(InputData(input1.to_dict()),
                            InputData(input2.to_dict()))
        input1 = input2.copy()
        if verbose and index.hour == 0:
            print(index.date(), f'SWE={model.snow_state.m_s:8.2f}')

    out = pd.DataFrame(model.output_list)
    return out


def swe_series(out):
    """Convenience: SWE indexed by datetime."""
    return out.set_index('date_time')['m_s']


def meltout(out, after_peak=True):
    """Datetime of melt-out, ignoring zero SWE before the seasonal peak."""
    swe = np.asarray(out['m_s'])
    start = int(np.argmax(swe)) if after_peak else 0
    gone = np.where(swe[start:] <= 0)[0]
    if not len(gone):
        return None
    return out['date_time'].iloc[start + gone[0]]


def state_at(out, date, columns=('m_s', 'z_s', 'rho')):
    """Model state at a given date. Returns a dict."""
    row = out.set_index('date_time').loc[pd.Timestamp(date)]
    return {c: float(row[c]) for c in columns}
