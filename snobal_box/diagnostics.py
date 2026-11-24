"""Comparison against reference runs, and SWE/depth identifiability."""

import numpy as np
import pandas as pd

__all__ = ['stats', 'find_offset', 'sensitivity_table', 'error_grid',
           'plot_error_grid', 'COMPARE_PAIRS']

def find_offset(ref, out, variable=('R_n', 'R_n_bar'), search=800):
    """Find the constant index offset between a reference .model file and a
    model run, for when the reference was produced from a different slice of
    the same forcing. Returns (offset, correlation)."""
    rcol, mcol = variable
    a = np.asarray(ref[rcol])
    b = np.asarray(out[mcol])
    best = (0, -np.inf)
    for lag in range(0, min(search, len(b) - 10)):
        n = min(len(a), len(b) - lag)
        if n < 100:
            break
        r = np.corrcoef(a[:n], b[lag:lag + n])[0, 1]
        if r > best[1]:
            best = (lag, r)
    return best


# reference .model column -> model output column
COMPARE_PAIRS = [
    ('m_s', 'm_s'), ('z_s', 'z_s'), ('rho', 'rho'), ('h2o', 'h2o'),
    ('R_n', 'R_n_bar'), ('H', 'H_bar'), ('L_v_E', 'L_v_E_bar'),
    ('delta_Q', 'delta_Q_bar'),
]


def stats(ref, out, pairs=None):
    """Agreement between a reference .model frame and a model run."""
    pairs = pairs or COMPARE_PAIRS
    n = min(len(ref), len(out))
    valid = ref['T_s_0'].values[:n] > -70   # -75 is the nodata sentinel
    rows = []
    for rcol, mcol in pairs:
        a = ref[rcol].values[:n][valid]
        b = out[mcol].values[:n][valid]
        rows.append({
            'variable': rcol,
            'ref_mean': a.mean(),
            'model_mean': b.mean(),
            'bias': (b - a).mean(),
            'rmse': np.sqrt(np.mean((b - a) ** 2)),
            'r': np.corrcoef(a, b)[0, 1],
        })
    return pd.DataFrame(rows).set_index('variable')


# ------------------------------------------------- SWE / depth identifiability

def state_at(out, date, columns=('m_s', 'z_s', 'rho')):
    """Model state at a given date. Returns a dict."""
    row = out.set_index('date_time').loc[pd.Timestamp(date)]
    return {c: float(row[c]) for c in columns}


def sensitivity_table(site, forcing, date, knobs=None, quiet=True):
    """Percent change in SWE, depth and bulk density at `date` for each knob.

    The point of this table is the RATIO of the SWE column to the depth
    column. Knobs that move both together are mass errors and leave bulk
    density alone. Knobs that move depth without SWE are density errors.
    """
    import warnings as _w
    knobs = knobs or [
        ('precip_scale', 1.10),
        ('threshold_shift', 1.0),
        ('rho_snow_scale', 1.20),
        ('T_offset', 1.0),
        ('solar_scale', 1.10),
        ('wind_scale', 1.30),
    ]
    with _w.catch_warnings():
        if quiet:
            _w.simplefilter('ignore')
        base = state_at(site.run(forcing), date)
        rows = []
        for knob, value in knobs:
            st = state_at(site.run(forcing, **{knob: value}), date)
            rows.append({
                'knob': f'{knob}={value}',
                'dSWE_%': 100 * (st['m_s'] - base['m_s']) / base['m_s'],
                'ddepth_%': 100 * (st['z_s'] - base['z_s']) / base['z_s'],
                'drho_%': 100 * (st['rho'] - base['rho']) / base['rho'],
            })
    df = pd.DataFrame(rows).set_index('knob')
    df['depth/SWE'] = df['ddepth_%'] / df['dSWE_%']
    return df


def error_grid(site, forcing, date, precip_scales=None, rho_scales=None,
               quiet=True, progress=True):
    """Map (precip_scale, rho_snow_scale) onto (SWE error, depth error).

    A mass knob and a density knob spanning a non-degenerate space means any
    observed pair of SWE and depth biases decomposes uniquely into the two.
    Returns a tidy DataFrame.
    """
    import warnings as _w
    precip_scales = precip_scales if precip_scales is not None \
        else np.arange(0.85, 1.16, 0.05)
    rho_scales = rho_scales if rho_scales is not None \
        else np.arange(0.60, 1.45, 0.10)

    rows = []
    with _w.catch_warnings():
        if quiet:
            _w.simplefilter('ignore')
        base = state_at(site.run(forcing), date)
        total = len(precip_scales) * len(rho_scales)
        for i, p in enumerate(precip_scales):
            for j, r in enumerate(rho_scales):
                st = state_at(site.run(forcing, precip_scale=p,
                                       rho_snow_scale=r), date)
                rows.append({
                    'precip_scale': p, 'rho_snow_scale': r,
                    'swe_err_%': 100 * (st['m_s'] - base['m_s']) / base['m_s'],
                    'depth_err_%': 100 * (st['z_s'] - base['z_s']) / base['z_s'],
                    'rho_err_%': 100 * (st['rho'] - base['rho']) / base['rho'],
                })
                if progress and len(rows) % 10 == 0:
                    print(f'  {len(rows)}/{total}', end='\r')
    if progress:
        print(' ' * 20, end='\r')
    return pd.DataFrame(rows), base


def plot_error_grid(grid, ax=None, levels_swe=None, levels_depth=None):
    """Contour SWE error and depth error over the (mass, density) plane."""
    import matplotlib.pyplot as plt
    if ax is None:
        _, ax = plt.subplots(figsize=(7.5, 6))

    piv_s = grid.pivot(index='rho_snow_scale', columns='precip_scale',
                       values='swe_err_%')
    piv_d = grid.pivot(index='rho_snow_scale', columns='precip_scale',
                       values='depth_err_%')
    X, Y = np.meshgrid(piv_s.columns.values, piv_s.index.values)

    cs = ax.contour(X, Y, piv_s.values, levels=levels_swe,
                    colors='#1f77b4', linewidths=1.2)
    ax.clabel(cs, fmt='%1.0f%%', fontsize=8)
    cd = ax.contour(X, Y, piv_d.values, levels=levels_depth,
                    colors='#d62728', linewidths=1.2, linestyles='--')
    ax.clabel(cd, fmt='%1.0f%%', fontsize=8)

    ax.plot(1.0, 1.0, 'ko', ms=7)
    ax.annotate('truth', (1.0, 1.0), textcoords='offset points',
                xytext=(8, 6), fontsize=9)
    ax.set_xlabel('precipitation mass scale')
    ax.set_ylabel('new snow density scale')
    ax.set_title('solid blue = SWE error,  dashed red = depth error')
    return ax
