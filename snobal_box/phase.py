"""The temperature to rain/snow lookup table, and repairs to precip files.

snobal does not decide rain versus snow. It is told, per timestep, what
fraction of precipitation is snow and what density that snow has. Those come
from a lookup table applied upstream of the model. PHASE_TABLE below was
reverse-engineered from the Atwater WY2024 precip file and is a starting
point, not an authority: substitute your own.
"""

import numpy as np
import pandas as pd

__all__ = ['PHASE_TABLE', 'apply_phase_table', 'shifted_phase_table',
           'phase_table_frame', 'repair_precip']

# Bins are (T_low, T_high, percent_snow, rho_snow), T in degrees C.
PHASE_TABLE = [
    (-np.inf, -5.0, 1.00,  75),
    (-5.0,    -3.0, 1.00, 100),
    (-3.0,    -1.0, 1.00, 150),
    (-1.0,    -0.5, 1.00, 175),
    (-0.5,     0.0, 0.75, 200),
    ( 0.0,     0.5, 0.25, 250),
    ( 0.5,  np.inf, 0.00,   0),
]


def apply_phase_table(T_pp_C, table=None):
    """Map precip temperature (C) to (percent_snow, rho_snow)."""
    table = table or PHASE_TABLE
    for lo, hi, pct, rho in table:
        if lo <= T_pp_C < hi:
            return pct, rho
    return 0.0, 0


def repair_precip(df, boundary_fill=175, verbose=True):
    """Fill rho_snow == 0 records that still have percent_snow > 0."""
    bad = (df['m_pp'] > 0) & (df['percent_snow'] > 0) & (df['rho_snow'] <= 0)
    if bad.any() and verbose:
        print(f'repairing {bad.sum()} precip record(s) with '
              f'percent_snow > 0 and rho_snow == 0:')
        print(df.loc[bad, ['m_pp', 'percent_snow', 'rho_snow', 'T_pp']])
    df = df.copy()
    df.loc[bad, 'rho_snow'] = boundary_fill
    return df


def shifted_phase_table(threshold_shift=0.0, table=None):
    """Return the phase table with every temperature threshold moved by
    `threshold_shift` degrees C. Positive moves the rain/snow boundary to a
    WARMER temperature, so more precipitation falls as snow."""
    table = table or PHASE_TABLE
    return [(lo + threshold_shift, hi + threshold_shift, pct, rho)
            for lo, hi, pct, rho in table]


def phase_table_frame(threshold_shift=0.0, table=None):
    """The phase table as a DataFrame, for printing or plotting."""
    return pd.DataFrame(shifted_phase_table(threshold_shift, table),
                        columns=['T_low_C', 'T_high_C',
                                 'percent_snow', 'rho_snow'])
