"""Scenario perturbations.

Every function here takes a forcing DataFrame and returns a forcing DataFrame.
Nothing touches model internals. That is deliberate: the scenario knobs, an
uploaded dataset and a coupled run all go through the same code path.
"""

import warnings

import numpy as np

from snobal_box._pysnobal.core.constants import (
    DATA_TSTEP, NORMAL_TSTEP, MEDIUM_TSTEP, SMALL_TSTEP, FREEZE)
from snobal_box.phase import apply_phase_table, shifted_phase_table

__all__ = ['perturb', 'snow_fraction']

def perturb(forcing, T_offset=0.0, precip_scale=1.0, solar_scale=1.0,
            lw_offset=0.0, wind_scale=1.0, threshold_shift=0.0,
            rho_snow_scale=1.0, rebuild_phase=None, phase_table=None,
            phase_shift=None):
    """Apply scenario perturbations to a forcing DataFrame.

    Operates entirely OUTSIDE the model: takes a forcing frame, returns a
    forcing frame. Nothing here touches snobal internals.

    T_offset         degrees C added to air temperature and precip temperature
    precip_scale     multiplies precipitation mass
    solar_scale      multiplies net solar (a crude albedo proxy)
    lw_offset        W/m2 added to incoming longwave
    wind_scale       multiplies wind speed

    threshold_shift  degrees C added to every temperature threshold in the
                     rain/snow lookup table. POSITIVE moves the rain/snow
                     boundary to a warmer temperature, so MORE precipitation
                     falls as snow. Negative gives more rain.
                     With the default table, the all-snow/some-rain boundary
                     sits at -0.5 C and the all-rain boundary at +0.5 C;
                     threshold_shift=+1 moves those to +0.5 and +1.5.

    rho_snow_scale   multiplies new snow density. Applied whether or not the
                     phase table is rebuilt.

    rebuild_phase    None (default) recomputes percent_snow and rho_snow from
                     the table whenever T_offset or threshold_shift is
                     nonzero. True always recomputes. False never does, which
                     lets you warm the air WITHOUT changing precipitation
                     phase -- useful for separating the energy-balance
                     response from the phase response.

    phase_table      override the table entirely, same shape as PHASE_TABLE.
    """
    if phase_shift is not None:            # old name, opposite documentation
        warnings.warn('phase_shift is renamed to threshold_shift; note the '
                      'documented sign was wrong, positive means MORE SNOW',
                      DeprecationWarning)
        threshold_shift = phase_shift

    df = forcing.copy()
    df['T_a'] = df['T_a'] + T_offset
    df['T_pp'] = df['T_pp'] + T_offset
    df['S_n'] = df['S_n'] * solar_scale
    df['I_lw'] = df['I_lw'] + lw_offset
    df['u'] = df['u'] * wind_scale
    df['m_pp'] = df['m_pp'] * precip_scale

    if rebuild_phase is None:
        rebuild_phase = (T_offset != 0.0) or (threshold_shift != 0.0)

    if rebuild_phase:
        table = shifted_phase_table(threshold_shift, phase_table)
        T_C = df['T_pp'] - FREEZE
        pct, rho = zip(*[apply_phase_table(t, table) for t in T_C])
        df['percent_snow'] = pct
        df['rho_snow'] = np.asarray(rho, dtype=float)

    if rho_snow_scale != 1.0:
        df['rho_snow'] = df['rho_snow'] * rho_snow_scale

    return df


def snow_fraction(forcing):
    """Fraction of total precipitation mass falling as snow."""
    total = forcing['m_pp'].sum()
    return float((forcing['m_pp'] * forcing['percent_snow']).sum() / total) \
        if total else float('nan')
