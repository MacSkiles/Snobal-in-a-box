"""Run configuration: one YAML file describes one site-season.

Everything a run needs that is not in the forcing lives here, so the same
forcing can be rerun with a different start time, precipitation scheme or
scenario set without touching the data. Relative paths are resolved against
the config file's folder.

See ``examples/`` for complete files and ``snobal-box init`` for a template.
"""

import copy
import os

import pandas as pd
import yaml

__all__ = ['load_config', 'resolve_config', 'DEFAULTS']

DEFAULTS = {
    'name': None,
    'site': {
        'elevation_m': None,
        'elevation_ft': None,
        'heights': None,           # {z_u, z_t, z_0, z_g} in m
        'height_file': None,       # or a snobal inheight .input file
        'initial_snow': None,      # {depth_m, density, temp_C}
        'snow_file': None,         # or a snobal snow.properties .input file
    },
    'forcing': {
        'format': None,            # input | csv
        # input (snobal .input files)
        'data_file': None,
        'precip_file': None,
        'start': None,             # 'YYYY-MM-DD HH:MM' | 'auto'
        'water_year': None,        # shorthand: start Oct 1 00:00 of WY-1
        'align_with': None,        # {file, columns, timezone, guess, search_days}
        # csv
        'file': None,
        'columns': 'snobal',       # preset name or {snobal_name: csv_name}
        'time_column': 'time',
        'max_gap_hours': 3,
        # both
        'timezone': 'UTC',         # of start / naive CSV times
        'output_timezone': None,   # report times in this zone (default: timezone)
        'temps_in_C': True,
        'window': None,            # {start, end} crop, in `timezone`
        'screen': True,            # replace out-of-range values
    },
    'precip': {},                  # see snobal_box.precip
    'model': {'max_h2o_vol': 0.01, 'max_z_s_0': 0.25},
    'scenarios': {},               # {label: {knob: value, precip: {...}}}
    'observations': None,          # {file, time_column, depth_column, ...}
    'output': {'dir': 'output', 'figures': True, 'write_forcing': False},
}

FT_TO_M = 0.3048


def _merge(base, over):
    out = copy.deepcopy(base)
    for k, v in (over or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict) \
                and k not in ('scenarios',):
            out[k] = _merge(out[k], v)
        else:
            out[k] = v
    return out


def load_config(path_or_dict, base_dir=None):
    """Read a YAML path (or take a dict) and return a resolved config."""
    if isinstance(path_or_dict, (str, os.PathLike)):
        with open(path_or_dict) as fh:
            raw = yaml.safe_load(fh) or {}
        base_dir = base_dir or os.path.dirname(os.path.abspath(path_or_dict))
    else:
        raw = path_or_dict
        base_dir = base_dir or os.getcwd()
    return resolve_config(raw, base_dir)


def resolve_config(raw, base_dir):
    cfg = _merge(DEFAULTS, raw)

    def p(x):
        if x is None or os.path.isabs(x):
            return x
        return os.path.normpath(os.path.join(base_dir, x))

    s, f = cfg['site'], cfg['forcing']
    for k in ('height_file', 'snow_file'):
        if s[k]:
            s[k] = p(s[k])
    for k in ('data_file', 'precip_file', 'file'):
        if f[k]:
            f[k] = p(f[k])
    if isinstance(f['align_with'], str):
        f['align_with'] = {'file': f['align_with']}
    if f['align_with']:
        f['align_with']['file'] = p(f['align_with']['file'])
    if cfg['observations'] and cfg['observations'].get('file'):
        cfg['observations']['file'] = p(cfg['observations']['file'])
    cfg['output']['dir'] = p(cfg['output']['dir'])

    # ---- checks that should fail before anything is read
    errors = []
    if not cfg['name']:
        errors.append('name: give the run a name')
    if s['elevation_m'] is None and s['elevation_ft'] is None:
        errors.append('site.elevation_m (or elevation_ft) is required')
    if s['elevation_m'] is None and s['elevation_ft'] is not None:
        s['elevation_m'] = float(s['elevation_ft']) * FT_TO_M
    if not (s['heights'] or s['height_file']):
        errors.append('site.heights {z_u, z_t, z_0, z_g} or '
                      'site.height_file is required')

    fmt = f['format'] or ('input' if f['data_file'] else
                          'csv' if f['file'] else None)
    if fmt == 'ipw':               # older name for the same format
        fmt = 'input'
    f['format'] = fmt
    if fmt == 'input':
        if not (f['data_file'] and f['precip_file']):
            errors.append('forcing.data_file and forcing.precip_file are '
                          'required for format input')
        if f['start'] is None and f['water_year'] is not None:
            f['start'] = f'{int(f["water_year"]) - 1}-10-01 00:00'
        if f['start'] is None:
            errors.append(
                'forcing.start is required for format input: the .input '
                'files have no timestamps. Give a start time or a '
                'water_year')
        if f['start'] == 'auto' and not f['align_with']:
            errors.append("forcing.start 'auto' needs forcing.align_with")
    elif fmt == 'csv':
        if not f['file']:
            errors.append('forcing.file is required for format csv')
    else:
        errors.append("forcing.format must be 'input' or 'csv'")

    for label, sc in cfg['scenarios'].items():
        if not isinstance(sc, dict):
            errors.append(f'scenarios.{label} must be a mapping of knobs')

    if errors:
        raise ValueError('config problems:\n  - ' + '\n  - '.join(errors))
    return cfg


def dump_config(cfg, path, extra=None):
    """Write the resolved config, plus anything learned while running (the
    start time actually used, alignment correlation), next to the output."""
    out = copy.deepcopy(cfg)
    if extra:
        out['_resolved'] = extra

    def clean(x):
        if isinstance(x, dict):
            return {k: clean(v) for k, v in x.items()}
        if isinstance(x, (list, tuple)):
            return [clean(v) for v in x]
        if isinstance(x, (pd.Timestamp,)):
            return str(x)
        if hasattr(x, 'item'):
            return x.item()
        return x
    with open(path, 'w') as fh:
        yaml.safe_dump(clean(out), fh, sort_keys=False)
