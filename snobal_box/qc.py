"""Flag forcing problems before they turn into model mysteries."""

import numpy as np

__all__ = ['qc_report']

def qc_report(forcing):
    """Flag the things that will bite before you spend a day on them."""
    issues = []
    n = len(forcing)

    gaps = forcing.index.to_series().diff().dropna().value_counts()
    if len(gaps) > 1:
        issues.append(f'irregular timestep: {dict(gaps)}')

    if forcing.isna().any().any():
        bad = forcing.isna().sum()
        issues.append(f'NaNs present: {dict(bad[bad > 0])}')

    # supersaturation: e_a above saturation makes h_le complain.
    # Use snobal's own sati()/satw() so this matches what the model checks.
    from snobal_box._pysnobal.point.libsnobal import sati
    es = np.array([sati(t) for t in forcing['T_a']])
    over = int((forcing['e_a'].values > es).sum())
    if over:
        issues.append(f'{over} records ({100 * over / n:.1f}%) with vapor '
                      'pressure above saturation')

    bad_precip = ((forcing['m_pp'] > 0) & (forcing['percent_snow'] > 0)
                  & (forcing['rho_snow'] <= 0)).sum()
    if bad_precip:
        issues.append(f'{bad_precip} precip records with percent_snow > 0 '
                      'and rho_snow == 0')

    if (forcing['u'] <= 0).any():
        issues.append('non-positive wind speed present')

    ranges = {
        'S_n': (0, 1400), 'I_lw': (100, 450), 'T_a': (220, 320),
        'e_a': (1, 2000), 'u': (0, 40), 'T_g': (230, 300),
    }
    for col, (lo, hi) in ranges.items():
        out = ((forcing[col] < lo) | (forcing[col] > hi)).sum()
        if out:
            issues.append(f'{col}: {out} values outside [{lo}, {hi}]')

    if issues:
        print('QC flags:')
        for i in issues:
            print('  -', i)
    else:
        print('QC: nothing flagged')
    return issues
