"""Season metrics for a SWE or depth series, and agreement with observations.

The same ``season_metrics`` is applied to model output and to observations, so
a comparison is always like with like: peak, peak date, half-melt date,
melt-out date and the depletion rate between peak and melt-out.
"""

import numpy as np
import pandas as pd

__all__ = ['season_metrics', 'compare', 'mass_balance', 'read_observations']


def season_metrics(series, zero=0.0, prefix=''):
    """Metrics of one snow season from a SWE (mm) or depth series.

    peak            maximum value
    peak_date       when it happened
    half_melt_date  first time after the peak the value is <= half the peak
    melt_out_date   first time after the peak the value is <= ``zero``
    depletion_days  peak to melt-out
    depletion_rate  peak / depletion_days (units per day)
    """
    s = series.dropna()
    if s.empty or s.max() <= zero:
        return {f'{prefix}peak': 0.0}
    peak_t = s.idxmax()
    peak = float(s.max())
    after = s.loc[peak_t:]
    half = after[after <= peak / 2]
    gone = after[after <= zero]
    out = {
        f'{prefix}peak': round(peak, 3),
        f'{prefix}peak_date': peak_t,
        f'{prefix}half_melt_date': half.index[0] if len(half) else pd.NaT,
        f'{prefix}melt_out_date': gone.index[0] if len(gone) else pd.NaT,
    }
    if len(gone):
        days = (gone.index[0] - peak_t) / pd.Timedelta('1D')
        out[f'{prefix}depletion_days'] = round(days, 2)
        out[f'{prefix}depletion_rate_per_day'] = \
            round(peak / days, 3) if days > 0 else np.nan
    return out


def compare(model, observed, name):
    """Bias, RMSE, r and timing differences of a model series against an
    observed one (any time spacing; the model is sampled at observation
    times)."""
    obs = observed.dropna()
    obs = obs[(obs.index >= model.index[0]) & (obs.index <= model.index[-1])]
    if obs.empty:
        return {}
    mod = model.reindex(model.index.union(obs.index)).interpolate(
        'time').reindex(obs.index)
    err = mod - obs
    out = {f'{name}_n_obs': int(len(obs)),
           f'{name}_bias': round(float(err.mean()), 4),
           f'{name}_rmse': round(float(np.sqrt((err ** 2).mean())), 4),
           f'{name}_r': round(float(np.corrcoef(mod, obs)[0, 1]), 4)
           if len(obs) > 2 else np.nan}
    m = season_metrics(model)
    o = season_metrics(obs)
    if 'peak' in o and o['peak'] > 0:
        out[f'{name}_peak_obs'] = o['peak']
        out[f'{name}_peak_error'] = round(m['peak'] - o['peak'], 4)
        for key in ('peak_date', 'half_melt_date', 'melt_out_date'):
            a, b = m.get(key, pd.NaT), o.get(key, pd.NaT)
            out[f'{name}_{key}_obs'] = b
            out[f'{name}_{key}_error_days'] = \
                round((a - b) / pd.Timedelta('1D'), 2) \
                if pd.notna(a) and pd.notna(b) else np.nan
    return out


def mass_balance(out, forcing, initial_swe):
    """Where the precipitation went, in mm, and the closure residual.

    snobal passes rain on bare ground straight through as snow water input,
    so closure counts all precipitation:

        initial SWE + precip - sublimation/evaporation - SWI - final SWE

    The residual should be a few millimetres at most (sub-step rounding). Rain on snow and on
    bare ground are reported separately for phase analyses.
    """
    snow_before = out['m_s'].shift(1).fillna(initial_swe) > 0
    f = forcing.reindex(out.index)
    m = f['m_pp'].fillna(0)
    snow = m * f['percent_snow'].fillna(0)
    rain = m - snow
    into_pack = m.sum()
    swi = out['ro_pred_sum'].sum()
    evap = out['E_s_sum'].sum()          # negative = mass lost
    final = out['m_s'].iloc[-1]
    return {
        'snowfall_mm': round(float(snow.sum()), 1),
        'rain_on_snow_mm': round(float(rain[snow_before].sum()), 1),
        'rain_on_bare_ground_mm': round(float(rain[~snow_before].sum()), 1),
        'snow_water_input_mm': round(float(swi), 1),
        'sublimation_evaporation_mm': round(float(-evap), 1),
        'final_swe_mm': round(float(final), 1),
        'closure_residual_mm': round(float(
            initial_swe + into_pack + evap - swi - final), 2),
    }


def read_observations(cfg, tz):
    """Read an observation CSV described by the ``observations`` config
    block. Returns a dict of Series in model units (depth m, SWE mm)."""
    raw = pd.read_csv(cfg['file'])
    t = pd.to_datetime(raw[cfg.get('time_column', 'time')])
    t = t.dt.tz_localize(cfg.get('timezone', tz)) if t.dt.tz is None \
        else t.dt.tz_convert(tz)
    out = {}
    to_m = {'m': 1.0, 'cm': 0.01, 'mm': 0.001, 'in': 0.0254}
    to_mm = {'mm': 1.0, 'cm': 10.0, 'm': 1000.0, 'in': 25.4}
    if cfg.get('depth_column'):
        s = pd.to_numeric(raw[cfg['depth_column']], errors='coerce')
        out['depth'] = pd.Series(
            s.values * to_m[cfg.get('depth_units', 'm')], index=t).sort_index()
    if cfg.get('swe_column'):
        s = pd.to_numeric(raw[cfg['swe_column']], errors='coerce')
        out['swe'] = pd.Series(
            s.values * to_mm[cfg.get('swe_units', 'mm')], index=t).sort_index()
    return out
