"""Config-driven runs: load forcing, build scenarios, run, summarize, save.

    import snobal_box as sb
    result = sb.run_config('examples/atwater_wy22_insitu.yaml')
    result.summary            # one row per scenario
    result.outputs['baseline']

or from a shell: ``snobal-box run examples/atwater_wy22_insitu.yaml``
"""

import math
import os
import warnings
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from snobal_box import __version__
from snobal_box._pysnobal.core.constants import FREEZE
from snobal_box.config import load_config, dump_config
from snobal_box.forcing import (read_input_forcing, read_csv_forcing,
                                read_input_heights, read_input_snow, MET_COLS)
from snobal_box.precip import apply_precip_scheme, resolve_auto
from snobal_box.align import find_start
from snobal_box.metrics import (season_metrics, compare, mass_balance,
                                read_observations)
from snobal_box.model import run as _run_core

__all__ = ['load_forcing', 'scenario_forcing', 'run_forcing', 'run_config',
           'RunResult', 'RANGES', 'SCENARIO_KNOBS']

# Plausible ranges in the units a human reads (temperatures in C).
RANGES = {
    'S_n': (0.0, 1400.0), 'I_lw': (100.0, 500.0), 'T_a': (-45.0, 45.0),
    'e_a': (1.0, 4000.0), 'u': (0.0, 40.0), 'T_g': (-30.0, 30.0),
}

SCENARIO_KNOBS = {
    'T_offset': 0.0,        # C added to air and precipitation temperature
    'solar_scale': 1.0,     # multiplies net solar
    'lw_offset': 0.0,       # W/m2 added to incoming longwave
    'wind_scale': 1.0,      # multiplies wind
    'vapor_scale': 1.0,     # multiplies vapor pressure
    'precip_scale': 1.0,    # multiplies precipitation mass (on top of config)
    'rho_snow_scale': 1.0,  # multiplies new snow density
    'threshold_shift': 0.0, # C; positive = rain/snow boundary warmer = MORE snow
    'rebuild_phase': None,  # False keeps the baseline phase under T_offset
    'precip': None,         # override any precip setting, e.g. {phase: dai2008}
}


# ---------------------------------------------------------------- loading

def screen(forcing, ranges=None, verbose=True):
    """Replace out-of-range meteorology with NaN, interpolate, report."""
    ranges = ranges or RANGES
    out, flagged = forcing.copy(), {}
    for col, (lo, hi) in ranges.items():
        if col not in out:
            continue
        human = out[col] - (FREEZE if col in ('T_a', 'T_g') else 0.0)
        bad = (human < lo) | (human > hi) | human.isna()
        if bad.any():
            flagged[col] = (int(bad.sum()),
                            sorted(human[bad].dropna().round(1).tolist())[:4])
            out.loc[bad, col] = np.nan
    if flagged:
        out[MET_COLS] = out[MET_COLS].interpolate(limit_direction='both')
        if verbose:
            for col, (n, ex) in flagged.items():
                print(f'  [QC] {col}: replaced {n} value(s), e.g. {ex}')
    return out, flagged


def _start_for_input(cfg, verbose):
    f = cfg['forcing']
    tz = f['timezone']
    info = {}
    if f['align_with']:
        a = f['align_with']
        ref = read_csv_forcing(a['file'], columns=a.get('columns', 'snobal'),
                               time_column=a.get('time_column', 'time'),
                               timezone=a.get('timezone', 'UTC'),
                               temps_in_C=a.get('temps_in_C', True))
        data = pd.read_csv(f['data_file'], sep=r'\s+', header=None,
                           usecols=[2]).iloc[:, 0].to_numpy()
        guess = a.get('guess') or (None if f['start'] == 'auto'
                                   else f['start'])
        if guess is not None:
            guess = pd.Timestamp(guess)
            if guess.tzinfo is None:
                guess = guess.tz_localize(tz)
        best, r, n = find_start(data, ref['T_a'], guess=guess,
                                search_days=a.get('search_days', 45))
        best = best.tz_convert(tz)
        info.update(aligned_start=str(best), align_r=round(r, 3),
                    align_overlap_h=n)
        if f['start'] == 'auto':
            start = best
            if verbose:
                print(f'  start found by alignment: {best} (r = {r:.3f})')
        else:
            start = pd.Timestamp(f['start'])
            start = start.tz_localize(tz) if start.tzinfo is None else start
            off = (best - start) / pd.Timedelta('1h')
            info['align_offset_h'] = off
            if verbose:
                tag = 'agrees' if off == 0 else f'is {off:+.0f} h off'
                print(f'  configured start {start} {tag} with the '
                      f'alignment ({best}, r = {r:.3f})')
            if abs(off) >= 2:
                warnings.warn(f'{cfg["name"]}: configured start differs from '
                              f'the air-temperature alignment by {off:+.0f} h')
    else:
        start = pd.Timestamp(f['start'])
        start = start.tz_localize(tz) if start.tzinfo is None else start
    info['start'] = str(start)
    return start, info


def load_forcing(cfg, verbose=True):
    """Read and screen the forcing for a resolved config. Precipitation
    phase/density are NOT applied yet. Returns (forcing, info)."""
    f = cfg['forcing']
    info = {}
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        if f['format'] == 'input':
            start, info = _start_for_input(cfg, verbose)
            forcing = read_input_forcing(f['data_file'], f['precip_file'],
                                       start, temps_in_C=f['temps_in_C'])
        else:
            forcing = read_csv_forcing(
                f['file'], columns=f['columns'],
                time_column=f['time_column'], timezone=f['timezone'],
                temps_in_C=f['temps_in_C'], max_gap_hours=f['max_gap_hours'])
    for w in caught:
        if verbose:
            print(f'  [read] {w.message}')
        info.setdefault('read_warnings', []).append(str(w.message))

    tz = f['output_timezone'] or f['timezone']
    forcing.index = forcing.index.tz_convert(tz)
    if f['window']:
        w = f['window']
        lo = pd.Timestamp(w['start']).tz_localize(tz) \
            if w.get('start') else forcing.index[0]
        hi = pd.Timestamp(w['end']).tz_localize(tz) \
            if w.get('end') else forcing.index[-1]
        forcing = forcing.loc[lo:hi]
    if f['screen']:
        forcing, flagged = screen(forcing, verbose=verbose)
        if flagged:
            info['screened'] = {k: v[0] for k, v in flagged.items()}
    if verbose:
        print(f'  {cfg["name"]}: {len(forcing)} hourly records, '
              f'{forcing.index[0]:%Y-%m-%d %H:%M %Z} to '
              f'{forcing.index[-1]:%Y-%m-%d %H:%M %Z}')
    return forcing, info


def site_state(cfg):
    """Measurement heights and initial snow from the config."""
    s = cfg['site']
    heights = dict(s['heights']) if s['heights'] \
        else read_input_heights(s['height_file'])
    for k in ('z_u', 'z_t', 'z_0', 'z_g'):
        if k not in heights:
            raise ValueError(f'site.heights is missing {k}')
    init = {'depth_m': 0.0, 'density': 300.0, 'temp_C': -4.0,
            'h2o_sat': 0.0}
    if s['snow_file']:
        init.update(read_input_snow(s['snow_file']))
    if s['initial_snow']:
        init.update(s['initial_snow'])
    t = init['temp_C'] + FREEZE
    t0 = init.get('temp_surface_C', init['temp_C']) + FREEZE
    snow = {'z_s': float(init['depth_m'] or 0.0),
            'rho': float(init['density']),
            't_s_0': t0, 't_s': t, 't_s_l': t,
            'h2o_sat': float(init['h2o_sat']), 'z_0': heights['z_0']}
    return heights, snow


# -------------------------------------------------------------- scenarios

def scenario_forcing(raw, precip_cfg=None, **knobs):
    """Forcing for one scenario: perturb the meteorology, then apply the
    precipitation scheme. Returns (forcing, notes)."""
    unknown = set(knobs) - set(SCENARIO_KNOBS)
    if unknown:
        raise ValueError(f'unknown scenario knob(s) {sorted(unknown)}; '
                         f'known: {sorted(SCENARIO_KNOBS)}')
    k = dict(SCENARIO_KNOBS, **knobs)
    precip = dict(precip_cfg or {})
    if k['precip']:
        precip.update(k['precip'])
    precip = resolve_auto(raw, precip)
    if k['threshold_shift']:
        precip['threshold_shift'] = precip.get('threshold_shift', 0.0) \
            + k['threshold_shift']
    notes = []

    df = raw.copy()
    df['T_a'] += k['T_offset']
    if 'T_pp' in df:
        df['T_pp'] += k['T_offset']
    df['S_n'] *= k['solar_scale']
    df['I_lw'] += k['lw_offset']
    df['u'] *= k['wind_scale']
    df['e_a'] *= k['vapor_scale']

    rebuild = k['rebuild_phase']
    if rebuild is None:
        rebuild = bool(k['T_offset'])
    needs_scheme = rebuild or precip['threshold_shift']
    if needs_scheme and precip['phase'] == 'file':
        precip['phase'] = 'table'
        notes.append("phase 'file' cannot respond to warming or a threshold "
                     "shift; used the lookup table for this scenario")
        if precip['density'] == 'file':
            precip['density'] = 'table'

    if rebuild or not k['T_offset']:
        df = apply_precip_scheme(df, precip)
    else:
        # warm the air but keep the baseline's rain/snow split and density:
        # isolates the energy-balance response from the phase response
        base = apply_precip_scheme(raw, precip)
        df['percent_snow'] = base['percent_snow']
        df['rho_snow'] = base['rho_snow']
        df['m_pp'] = base['m_pp']
        df['T_pp'] = base['T_pp'] + k['T_offset']
        notes.append('phase and density held at baseline (rebuild_phase: '
                     'false)')

    df['m_pp'] *= k['precip_scale']
    df['rho_snow'] *= k['rho_snow_scale']
    notes.append(f"phase={precip['phase']} density={precip['density']} "
                 f"on {precip['phase_temperature']}")
    return df, notes


# -------------------------------------------------------------------- run

_FALLBACKS = []
_PATCHED = False


def _patch_hle1():
    """Substitute neutral-stability fluxes when snobal's turbulent-transfer
    iteration fails to converge, instead of killing the run. Each use is
    counted and reported. A handful in a season is noise; hundreds mean the
    forcing needs a look."""
    global _PATCHED
    if _PATCHED:
        return
    from snobal_box._pysnobal.point import libsnobal as L
    original = L.hle1

    def hle1_safe(press, ta, ts, za, ea, es, zq, u, zu, z0):
        h, le, e, status = original(press, ta, ts, za, ea, es, zq, u, zu, z0)
        if status == 0:
            return h, le, e, status
        es = min(es, L.sati(ts))
        ea = min(ea, L.satw(ta))
        d0 = 2 * L.PAESCHKE * z0 / 3
        ltsm, ltsh = math.log((zu - d0) / z0), math.log((za - d0) / z0)
        ltsv = math.log((zq - d0) / z0)
        ta_pot = ta + L.DALR * za
        dens = L.GAS_DEN(press, L.MOL_AIR,
                         L.VIR_TEMP(math.sqrt(ta_pot * ts),
                                    math.sqrt(ea * es), press))
        factor = L.VON_KARMAN * (L.VON_KARMAN * u / ltsm) * dens
        e_n = (L.spec_hum(ea, press) - L.spec_hum(es, press)) * factor \
            * L.AV / ltsv
        h_n = (ta_pot - ts) * factor * L.CP_AIR * L.AH / ltsh
        xlh = L.LH_VAP(ts) + (L.LH_FUS(ts) if ts <= L.FREEZE else 0.0)
        _FALLBACKS.append(1)
        return h_n, xlh * e_n, e_n, 0

    L.hle1 = hle1_safe
    _PATCHED = True


def run_forcing(forcing, heights, snow, elevation_m, max_h2o_vol=0.01,
                max_z_s_0=0.25):
    """Run snobal over a ready forcing frame. Returns (output, n_fallbacks)."""
    _patch_hle1()
    _FALLBACKS.clear()
    cols = ['S_n', 'I_lw', 'T_a', 'e_a', 'u', 'T_g', 'm_pp',
            'percent_snow', 'rho_snow', 'T_pp']
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        out = _run_core(forcing[cols], snow, heights, elevation=elevation_m,
                        max_h2o_vol=max_h2o_vol, max_z_s_0=max_z_s_0)
    out = out.set_index('date_time')
    out.index.name = 'time'
    return out, len(_FALLBACKS)


@dataclass
class RunResult:
    config: dict
    summary: pd.DataFrame
    outputs: dict = field(default_factory=dict)
    forcings: dict = field(default_factory=dict)
    observations: dict = field(default_factory=dict)
    info: dict = field(default_factory=dict)
    paths: dict = field(default_factory=dict)


def _summarize(label, out, forcing, snow, n_fb, notes, obs, heights):
    init_swe = snow['z_s'] * snow['rho']
    row = {'scenario': label}
    m = season_metrics(out['m_s'], prefix='swe_')
    d = season_metrics(out['z_s'], prefix='depth_')
    row.update(m)
    row['depth_peak_m'] = d.get('depth_peak')
    row['depth_peak_date'] = d.get('depth_peak_date')
    row['depth_melt_out_date'] = d.get('depth_melt_out_date')
    row.update(mass_balance(out, forcing, init_swe))
    m_pp = forcing['m_pp'].sum()
    row['precip_mm'] = round(float(m_pp), 1)
    row['snow_fraction'] = round(
        float((forcing['m_pp'] * forcing['percent_snow']).sum() / m_pp), 4) \
        if m_pp else np.nan
    row['melt_mm'] = round(float(out['melt_sum'].sum()), 1)
    row['neutral_stability_steps'] = n_fb
    row['max_depth_over_wind_height'] = round(
        float(out['z_s'].max() / heights['z_u']), 3)
    for name, series in obs.items():
        col = 'z_s' if name == 'depth' else 'm_s'
        row.update(compare(out[col], series, name))
    row['notes'] = '; '.join(notes)
    return row


def run_config(config, verbose=True, write=True, only=None):
    """Run the baseline and every scenario in a config. ``config`` is a YAML
    path or a dict. ``only`` limits to a list of scenario labels."""
    cfg = load_config(config)
    if verbose:
        print(f'== {cfg["name"]}  (snobal-box {__version__})')
    raw, info = load_forcing(cfg, verbose=verbose)
    heights, snow = site_state(cfg)
    obs = read_observations(cfg['observations'],
                            cfg['forcing']['output_timezone']
                            or cfg['forcing']['timezone']) \
        if cfg['observations'] else {}

    scenarios = {'baseline': {}}
    scenarios.update(cfg['scenarios'] or {})
    if only:
        scenarios = {k: v for k, v in scenarios.items() if k in only}

    rows, outputs, forcings = [], {}, {}
    for label, knobs in scenarios.items():
        forcing, notes = scenario_forcing(raw, cfg['precip'], **knobs)
        out, n_fb = run_forcing(forcing, heights, snow,
                                cfg['site']['elevation_m'],
                                **cfg['model'])
        row = _summarize(label, out, forcing, snow, n_fb, notes, obs,
                         heights)
        rows.append(row)
        outputs[label], forcings[label] = out, forcing
        if verbose:
            print(f'  {label:>24}: peak SWE {row["swe_peak"]:7.1f} mm on '
                  f'{row["swe_peak_date"]:%Y-%m-%d}, melt-out '
                  f'{_fmt(row.get("swe_melt_out_date"))}, snow fraction '
                  f'{row["snow_fraction"]:.2f}'
                  + (f', {n_fb} neutral-stability steps' if n_fb else ''))
            if row['max_depth_over_wind_height'] > 0.8:
                print(f'  {"":>24}  [!] depth reaches '
                      f'{row["max_depth_over_wind_height"]:.0%} of the wind '
                      'sensor height')

    summary = pd.DataFrame(rows).set_index('scenario')
    summary.insert(0, 'run', cfg['name'])
    result = RunResult(cfg, summary, outputs, forcings, obs, info)
    if write:
        result.paths = save(result, verbose=verbose)
    return result


def _fmt(t):
    return f'{t:%Y-%m-%d}' if isinstance(t, pd.Timestamp) and pd.notna(t) \
        else 'none'


def _slug(s):
    s = str(s).replace('+', 'plus').replace('-', 'minus')
    s = ''.join(c if c.isalnum() or c == '.' else '_' for c in s)
    while '__' in s:
        s = s.replace('__', '_')
    return s.strip('_').lower()


def save(result, verbose=True):
    cfg = result.config
    outdir = os.path.join(cfg['output']['dir'], _slug(cfg['name']))
    os.makedirs(outdir, exist_ok=True)
    paths = {}
    p = os.path.join(outdir, 'summary.csv')
    result.summary.to_csv(p)
    paths['summary'] = p
    for label, out in result.outputs.items():
        p = os.path.join(outdir, f'output_{_slug(label)}.csv')
        out.to_csv(p, float_format='%.5g')
        paths[f'output_{label}'] = p
        if cfg['output']['write_forcing']:
            p = os.path.join(outdir, f'forcing_{_slug(label)}.csv')
            result.forcings[label].to_csv(p, float_format='%.5g')
    p = os.path.join(outdir, 'config_resolved.yaml')
    dump_config(cfg, p, dict(result.info, snobal_box_version=__version__))
    paths['config'] = p
    if cfg['output']['figures']:
        from snobal_box import plots
        paths.update(plots.save_figures(result, outdir))
    if verbose:
        print(f'  wrote {len(paths)} file(s) to {outdir}')
    return paths
