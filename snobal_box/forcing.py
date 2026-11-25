"""Read forcing from either supported format into one hourly DataFrame.

Two formats:

* ``input``: the snobal ``.input`` pair (data file + precip file). These
  carry an hourly counter, not timestamps, so the start time comes from the
  config, or is found by matching air temperature against a timestamped
  reference (see :mod:`snobal_box.align`).
* ``csv``: one file with a time column. This is the format for uploads. Column
  names are mapped through a preset or an explicit ``columns`` dict.

Whatever the source, the result has the same columns, in snobal's units
(temperatures in kelvin) and a timezone-aware hourly index:

    S_n, I_lw, T_a, e_a, u, T_g, m_pp, percent_snow, rho_snow, T_pp

``percent_snow``, ``rho_snow`` and ``T_pp`` may be absent from a CSV. They
are then filled by the precipitation scheme (:mod:`snobal_box.precip`).
"""

import re
import warnings

import numpy as np
import pandas as pd

from snobal_box._pysnobal.core.constants import FREEZE

__all__ = ['read_input_forcing', 'read_csv_forcing', 'COLUMN_PRESETS',
           'MET_COLS', 'PHASE_COLS', 'FORCING_COLS', 'TEMP_COLS']

MET_COLS = ['S_n', 'I_lw', 'T_a', 'e_a', 'u', 'T_g']
PHASE_COLS = ['percent_snow', 'rho_snow', 'T_pp']
FORCING_COLS = MET_COLS + ['m_pp'] + PHASE_COLS
TEMP_COLS = ['T_a', 'T_g', 'T_pp']

# Map from snobal names to the column names used in a given CSV layout.
COLUMN_PRESETS = {
    # snobal's own names; the format of the upload template
    'snobal': {c: c for c in FORCING_COLS},
    # HRRR point extractions (jordan_bilinear_forcing_*.csv)
    'hrrr_extract': {
        'S_n': 'net_solar_Wm-2',
        'I_lw': 'downwelling_thermal_Wm-2',
        'T_a': 'temp_air_degC',
        'e_a': 'vapor_pressure_Pa',
        'u': 'wind_speed_ms-1',
        'T_g': 'temp_ground_degC',
        'm_pp': 'precip_mass_mm',
        'T_pp': 'precip_temp_degC',
        'percent_snow': 'snow_precip_fraction',
        'rho_snow': 'snow_precip_density_kgm-3',
    },
}


# ------------------------------------------------------------ .input files

def _read_fixed_fields(path, names, fill_missing=None):
    """Read whitespace-separated records that should have len(names) fields.

    Some exported files leave a field empty, which a plain whitespace split
    silently turns into a column shift (Atwater WY2024: the first 24 precip
    rows have no precipitation mass, so density received the temperature).
    When a row is one field short and splitting on single spaces gives the
    right count, the empty field is read as missing. Anything else is an
    error naming the line.
    """
    rows, repaired = [], []
    with open(path) as fh:
        for i, line in enumerate(fh, start=1):
            line = line.rstrip('\r\n')
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) != len(names):
                single = line.strip().split(' ')
                if len(single) == len(names) and single.count('') >= 1:
                    parts = [p if p != '' else 'nan' for p in single]
                    repaired.append(i)
                else:
                    raise ValueError(
                        f'{path}: line {i} has {len(parts)} fields, expected '
                        f'{len(names)}: {line!r}')
            rows.append([float(p) for p in parts])
    df = pd.DataFrame(rows, columns=names)
    if repaired:
        cols = [c for c in names if df.loc[[r - 1 for r in repaired], c]
                .isna().any()] if len(df) else []
        msg = (f'{path}: {len(repaired)} line(s) had an empty field '
               f'({", ".join(cols)}), first at line {repaired[0]}')
        if fill_missing:
            for c in cols:
                if c in fill_missing:
                    df[c] = df[c].fillna(fill_missing[c])
            msg += (f'; filled with {fill_missing}. Empty rows can mean '
                    'the column is shifted in time, so check the file')
        warnings.warn(msg)
    return df


def read_input_forcing(data_file, precip_file, start, temps_in_C=True,
                     tstep_minutes=60):
    """Read a snobal .input pair. ``start`` must be a tz-aware Timestamp."""
    data = _read_fixed_fields(data_file, MET_COLS)
    precip = _read_fixed_fields(
        precip_file, ['time', 'm_pp', 'percent_snow', 'rho_snow', 'T_pp'],
        fill_missing={'m_pp': 0.0}).drop(columns='time')
    if len(data) != len(precip):
        raise ValueError(f'{data_file} has {len(data)} records but '
                         f'{precip_file} has {len(precip)}')
    df = pd.concat([data, precip], axis=1)
    df.index = pd.date_range(start=start, periods=len(df),
                             freq=f'{tstep_minutes}min', name='time')
    if temps_in_C:
        for c in TEMP_COLS:
            df[c] = df[c] + FREEZE
    return df


def read_input_heights(path):
    """time z_u z_t z_0 z_g"""
    rec = _read_fixed_fields(path, ['t', 'z_u', 'z_t', 'z_0', 'z_g']).iloc[0]
    return {k: float(rec[k]) for k in ('z_u', 'z_t', 'z_0', 'z_g')}


def read_input_snow(path):
    """time z_s rho T_s_0 T_s h2o_sat. Temperatures returned in C."""
    rec = _read_fixed_fields(
        path, ['t', 'z_s', 'rho', 'T_s_0', 'T_s', 'h2o_sat']).iloc[0]
    return {'depth_m': float(rec['z_s']), 'density': float(rec['rho']),
            'temp_surface_C': float(rec['T_s_0']),
            'temp_C': float(rec['T_s']), 'h2o_sat': float(rec['h2o_sat'])}


# --------------------------------------------------------------------- CSV

def read_csv_forcing(path, columns='snobal', time_column='time',
                     timezone='UTC', temps_in_C=True, max_gap_hours=3):
    """Read a timestamped CSV.

    ``columns`` is a preset name from COLUMN_PRESETS or a dict mapping snobal
    names to the CSV's names. ``timezone`` applies to naive timestamps; times
    that already carry an offset are converted to it.

    Short gaps (up to ``max_gap_hours``) in the hourly record are filled:
    meteorology by linear interpolation, precipitation with zero. Longer gaps
    are an error, because snobal has no way to skip time.
    """
    mapping = COLUMN_PRESETS[columns] if isinstance(columns, str) \
        else dict(columns)
    raw = pd.read_csv(path)
    if time_column not in raw:
        raise ValueError(f'{path}: no time column {time_column!r}; columns '
                         f'are {list(raw.columns)}')

    required = MET_COLS + ['m_pp']
    missing = [f'{k} (looked for {mapping.get(k, k)!r})' for k in required
               if mapping.get(k, k) not in raw]
    if missing:
        raise ValueError(f'{path}: missing required column(s): '
                         + ', '.join(missing))

    t = pd.to_datetime(raw[time_column])
    t = t.dt.tz_localize(timezone) if t.dt.tz is None \
        else t.dt.tz_convert(timezone)

    df = pd.DataFrame(index=pd.DatetimeIndex(t, name='time'))
    for k in FORCING_COLS:
        src = mapping.get(k, k)
        if src in raw:
            df[k] = pd.to_numeric(raw[src], errors='coerce').values

    if df.index.duplicated().any():
        raise ValueError(f'{path}: duplicated timestamps, e.g. '
                         f'{df.index[df.index.duplicated()][0]}')
    df = df.sort_index()
    full = pd.date_range(df.index[0], df.index[-1], freq='h', name='time')
    if len(full) != len(df):
        gaps = full.difference(df.index)
        runs = (gaps.to_series().diff() != pd.Timedelta('1h')).cumsum()
        longest = int(runs.value_counts().max())
        if longest > max_gap_hours:
            raise ValueError(
                f'{path}: {len(gaps)} missing hour(s), longest run '
                f'{longest} h (max_gap_hours={max_gap_hours}); first missing '
                f'{gaps[0]}')
        warnings.warn(f'{path}: filled {len(gaps)} missing hour(s)')
        df = df.reindex(full)
        df['m_pp'] = df['m_pp'].fillna(0.0)

    if temps_in_C:
        for c in TEMP_COLS:
            if c in df:
                df[c] = df[c] + FREEZE
    met = [c for c in MET_COLS + ['T_pp'] if c in df]
    df[met] = df[met].interpolate(limit=max_gap_hours,
                                  limit_direction='both')
    for c in ('percent_snow', 'rho_snow'):
        if c in df:
            df[c] = df[c].fillna(0.0)
    return df
