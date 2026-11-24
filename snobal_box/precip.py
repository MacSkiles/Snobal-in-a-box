"""Precipitation phase and new snow density schemes.

snobal never decides rain versus snow. Every hour it is told how much water
fell (``m_pp``), what fraction of it was snow (``percent_snow``), how dense
that new snow is (``rho_snow``) and how warm the precipitation was
(``T_pp``). This module produces the last three from the meteorology, so the
choice is explicit and swappable instead of baked into a file.

A scheme is configured as a dict, normally the ``precip`` block of a run
config::

    precip:
      scale: 1.0                 # multiplies m_pp
      phase: table               # file | table | threshold | linear | dai2008
      phase_params: {}           # scheme specific, see PHASE_SCHEMES
      phase_temperature: T_pp    # T_pp | T_a | dewpoint
      density: table             # file | table | constant | hedstrom_pomeroy
      density_params: {}
      threshold_shift: 0.0       # C; positive moves the boundary warmer = MORE snow

``file`` means "use the values already in the forcing" and is only possible
when the forcing carries them (IPW precip files and the HRRR extractions do).
``auto`` (the default) is ``file`` when the forcing has the column and
``table`` when it does not.

Adding a scheme is one function plus one line in PHASE_SCHEMES or
DENSITY_SCHEMES.
"""

import warnings

import numpy as np
import pandas as pd

from snobal_box._pysnobal.core.constants import FREEZE
from snobal_box.phase import PHASE_TABLE

__all__ = ['apply_precip_scheme', 'PHASE_SCHEMES', 'DENSITY_SCHEMES',
           'DEFAULT_PRECIP', 'dewpoint_C', 'precip_partition']

DEFAULT_PRECIP = dict(scale=1.0, phase='auto', phase_params={},
                      phase_temperature='T_pp', density='auto',
                      density_params={}, threshold_shift=0.0,
                      rain_temperature='keep')


def dewpoint_C(e_a_Pa):
    """Dew point (C) from vapor pressure (Pa); Magnus form over water
    (Alduchov and Eskridge 1996 coefficients)."""
    lg = np.log(np.maximum(np.asarray(e_a_Pa, dtype=float), 1e-3) / 610.94)
    return 243.04 * lg / (17.625 - lg)


# ------------------------------------------------------------ phase schemes
# Each takes precipitation temperature in C (numpy array) and returns the
# snow fraction in [0, 1].

def _phase_table(T, table=None):
    table = table or PHASE_TABLE
    out = np.zeros_like(T)
    for lo, hi, pct, _ in table:
        out[(T >= lo) & (T < hi)] = pct
    return out


def _phase_threshold(T, T_threshold=0.0):
    """All snow below the threshold, all rain at or above it."""
    return (T < T_threshold).astype(float)


def _phase_linear(T, T_snow=-1.0, T_rain=3.0):
    """All snow at or below T_snow, all rain at or above T_rain, linear
    in between."""
    if T_rain <= T_snow:
        raise ValueError('linear phase scheme needs T_rain > T_snow')
    return np.clip((T_rain - T) / (T_rain - T_snow), 0.0, 1.0)


def _phase_dai2008(T, a=-48.2292, b=0.7205, c=1.1662, d=1.0223):
    """Hyperbolic tangent fit of snow frequency to air temperature over land
    (Dai 2008, GRL, doi:10.1029/2008GL033295). Defaults are the annual
    global land coefficients; check them against the paper for your use."""
    return np.clip(a * (np.tanh(b * (T - c)) - d) / 100.0, 0.0, 1.0)


PHASE_SCHEMES = {
    'table': _phase_table,
    'threshold': _phase_threshold,
    'linear': _phase_linear,
    'dai2008': _phase_dai2008,
}


# ---------------------------------------------------------- density schemes
# Each takes precipitation temperature in C and returns new snow density in
# kg/m3 (used only where snow fraction > 0).

def _density_table(T, table=None):
    table = table or PHASE_TABLE
    out = np.zeros_like(T)
    for lo, hi, _, rho in table:
        out[(T >= lo) & (T < hi)] = rho
    # the table's warmest bin carries density 0 (all rain); a scheme that
    # still calls some of that snow needs a density, so use the warmest
    # nonzero one
    warmest = [rho for *_, rho in table if rho > 0][-1]
    return np.where(out > 0, out, warmest)


def _density_constant(T, rho=100.0):
    return np.full_like(T, float(rho))


def _density_hedstrom_pomeroy(T, max_T=0.0):
    """Hedstrom and Pomeroy (1998): 67.92 + 51.25 exp(T / 2.59). The fit is
    for sub-freezing snowfall, so T is capped at ``max_T`` (C)."""
    T = np.minimum(T, max_T)
    return 67.92 + 51.25 * np.exp(T / 2.59)


DENSITY_SCHEMES = {
    'table': _density_table,
    'constant': _density_constant,
    'hedstrom_pomeroy': _density_hedstrom_pomeroy,
}


# ------------------------------------------------------------------- apply

def _phase_temperature_C(df, which):
    if which == 'T_pp':
        if 'T_pp' not in df or df['T_pp'].isna().all():
            return (df['T_a'] - FREEZE).to_numpy()
        return (df['T_pp'] - FREEZE).to_numpy()
    if which == 'T_a':
        return (df['T_a'] - FREEZE).to_numpy()
    if which == 'dewpoint':
        return dewpoint_C(df['e_a'].to_numpy())
    raise ValueError(f'phase_temperature must be T_pp, T_a or dewpoint, '
                     f'not {which!r}')


def resolve_auto(forcing, precip=None):
    """Fill defaults and turn 'auto' into 'file' or 'table'."""
    cfg = dict(DEFAULT_PRECIP, **(precip or {}))
    for key, col in (('phase', 'percent_snow'), ('density', 'rho_snow')):
        if cfg[key] == 'auto':
            has = col in forcing and not forcing[col].isna().all()
            cfg[key] = 'file' if has else 'table'
    return cfg


def apply_precip_scheme(forcing, precip=None, verbose=False):
    """Return a copy of ``forcing`` with m_pp scaled and percent_snow,
    rho_snow and T_pp set according to ``precip`` (see module docstring)."""
    cfg = resolve_auto(forcing, precip)
    df = forcing.copy()

    if 'T_pp' not in df or df['T_pp'].isna().all():
        df['T_pp'] = df['T_a']           # no precip temperature: use air
    if cfg['rain_temperature'] == 'T_a':
        df['T_pp'] = df['T_a']
    elif cfg['rain_temperature'] == 'dewpoint':
        df['T_pp'] = dewpoint_C(df['e_a'].to_numpy()) + FREEZE

    # evaluating the scheme at (T - shift) moves every threshold by +shift
    T = _phase_temperature_C(df, cfg['phase_temperature']) \
        - float(cfg['threshold_shift'])

    if cfg['phase'] == 'file':
        if 'percent_snow' not in df:
            raise ValueError("precip phase is 'file' but the forcing has no "
                             "percent_snow column; choose a phase scheme")
        if cfg['threshold_shift']:
            raise ValueError("threshold_shift needs a phase scheme, not "
                             "'file'")
    else:
        fn = PHASE_SCHEMES[cfg['phase']]
        df['percent_snow'] = fn(T, **cfg['phase_params'])

    if cfg['density'] == 'file':
        if 'rho_snow' not in df:
            raise ValueError("precip density is 'file' but the forcing has "
                             "no rho_snow column; choose a density scheme")
        bad = (df['m_pp'] > 0) & (df['percent_snow'] > 0) \
            & ~(df['rho_snow'] > 0)
        if bad.any():
            # a file density of 0 where the phase now says snow; snobal
            # would divide by it. Fill from the table and say so.
            fill = _density_table(T[bad.to_numpy()])
            df.loc[bad, 'rho_snow'] = fill
            if verbose:
                print(f'  [precip] {int(bad.sum())} snowy hour(s) had no '
                      'density in the file; filled from the lookup table')
    else:
        fn = DENSITY_SCHEMES[cfg['density']]
        df['rho_snow'] = fn(T, **cfg['density_params'])

    df['rho_snow'] = np.where(df['percent_snow'] > 0, df['rho_snow'], 0.0)
    df['m_pp'] = df['m_pp'] * float(cfg['scale'])
    return df


def precip_partition(forcing, swe=None):
    """Totals (mm) of precipitation, snowfall and rainfall. With a SWE
    series aligned to the forcing, also rain falling on snow."""
    m = forcing['m_pp']
    snow = m * forcing['percent_snow']
    rain = m - snow
    out = {'precip_mm': float(m.sum()), 'snowfall_mm': float(snow.sum()),
           'rainfall_mm': float(rain.sum()),
           'snow_fraction': float(snow.sum() / m.sum()) if m.sum() else
           float('nan')}
    if swe is not None:
        on = swe.reindex(forcing.index).shift(1).fillna(0) > 0
        out['rain_on_snow_mm'] = float(rain[on].sum())
    return out
