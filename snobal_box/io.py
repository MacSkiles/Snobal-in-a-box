"""Readers for classic IPW-style ASCII snobal files."""

import os
import warnings

import pandas as pd

from snobal_box._pysnobal.core.constants import (
    DATA_TSTEP, NORMAL_TSTEP, MEDIUM_TSTEP, SMALL_TSTEP, FREEZE)

__all__ = ['read_forcing', 'read_initial_snow', 'read_measurement_heights',
           'read_reference_output', 'DATA_COLS', 'PRECIP_COLS', 'OUT_COLS']

DATA_COLS = ['S_n', 'I_lw', 'T_a', 'e_a', 'u', 'T_g']
PRECIP_COLS = ['time', 'm_pp', 'percent_snow', 'rho_snow', 'T_pp']

# 26-column merged IPW snobal output (energy/mass + snow properties)
OUT_COLS = [
    'time', 'R_n', 'H', 'L_v_E', 'G', 'M', 'delta_Q', 'G_0', 'delta_Q_0',
    'cc_s_0', 'cc_s_l', 'cc_s', 'E_s', 'melt', 'ro_predict',
    'z_s_0', 'z_s_l', 'z_s', 'rho', 'm_s_0', 'm_s_l', 'm_s', 'h2o',
    'T_s_0', 'T_s_l', 'T_s']


def read_forcing(data_file, precip_file, start_date='2024-01-01',
                 tstep_minutes=60, temps_in_C=True):
    """Read the .input data and precip files into one hourly DataFrame."""
    data = pd.read_csv(data_file, sep=r'\s+', header=None, names=DATA_COLS)
    precip = pd.read_csv(precip_file, sep=r'\s+', header=None,
                         names=PRECIP_COLS).drop(columns='time')

    if len(data) != len(precip):
        raise ValueError(
            f'data has {len(data)} records, precip has {len(precip)}')

    df = pd.concat([data, precip], axis=1)
    df.index = pd.date_range(start=start_date, periods=len(df),
                             freq=f'{tstep_minutes}min')

    if temps_in_C:
        for col in ['T_a', 'T_g', 'T_pp']:
            df[col] = df[col] + FREEZE

    return df


def read_initial_snow(snow_file, z_0, temps_in_C=True):
    """Read the initial snowcover record."""
    rec = pd.read_csv(snow_file, sep=r'\s+', header=None).iloc[0]
    t_s_0, t_s = float(rec[3]), float(rec[4])
    if temps_in_C:
        t_s_0 += FREEZE
        t_s += FREEZE
    return {
        'z_s': float(rec[1]),
        'rho': float(rec[2]),
        't_s_0': t_s_0,
        't_s': t_s,
        't_s_l': t_s,
        'h2o_sat': float(rec[5]),
        'z_0': z_0,
    }


def read_measurement_heights(mh_file):
    """time, z_u, z_T, z_0, z_g"""
    rec = pd.read_csv(mh_file, sep=r'\s+', header=None).iloc[0]
    return {'z_u': float(rec[1]), 'z_t': float(rec[2]),
            'z_0': float(rec[3]), 'z_g': float(rec[4])}


def read_reference_output(model_file):
    """Read a reference .model file. Returns None if it is missing or is a
    degenerate stub (a failed C run often leaves a single all-nodata row)."""
    if model_file is None or not os.path.exists(model_file):
        return None
    ref = pd.read_csv(model_file, sep=r'\s+', header=None, names=OUT_COLS)
    if len(ref) < 2 or (ref['T_s_0'] <= -70).all():
        warnings.warn(
            f'{os.path.basename(model_file)} has {len(ref)} record(s) and no '
            'valid snow temperatures; treating as absent')
        return None
    return ref
