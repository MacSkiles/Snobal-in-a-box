"""Command line: ``snobal-box run|check|init``."""

import argparse
import os
import shutil
import sys
import warnings

import pandas as pd

__all__ = ['main']

TEMPLATES = os.path.join(os.path.dirname(__file__), 'templates')


def _run(args):
    from snobal_box.runner import run_config
    tables = []
    for path in args.config:
        res = run_config(path, verbose=not args.quiet,
                         only=args.scenario or None)
        tables.append(res.summary)
        print()
    if args.combined and tables:
        pd.concat(tables).to_csv(args.combined)
        print(f'combined summary: {args.combined}')


def _check(args):
    from snobal_box.config import load_config
    from snobal_box.runner import load_forcing, site_state, scenario_forcing
    from snobal_box.qc import qc_report
    for path in args.config:
        cfg = load_config(path)
        print(f'== {cfg["name"]}')
        raw, info = load_forcing(cfg)
        heights, snow = site_state(cfg)
        print(f'  heights {heights}')
        print(f'  initial snow depth {snow["z_s"]} m, density {snow["rho"]}')
        forcing, notes = scenario_forcing(raw, cfg['precip'])
        print(f'  precip: {"; ".join(notes)}')
        m = forcing['m_pp']
        print(f'  precip {m.sum():.0f} mm, '
              f'{(m * forcing["percent_snow"]).sum() / max(m.sum(), 1e-9):.0%}'
              ' as snow')
        qc_report(forcing)
        print()


def _init(args):
    os.makedirs(args.dir, exist_ok=True)
    for name in ('my_site.yaml', 'forcing_template.csv'):
        dst = os.path.join(args.dir, name)
        if os.path.exists(dst) and not args.force:
            print(f'exists, not overwritten: {dst}')
            continue
        shutil.copy(os.path.join(TEMPLATES, name), dst)
        print(f'wrote {dst}')
    print('\nFill in my_site.yaml, put your hourly forcing in the CSV '
          '(one row per hour, column meanings are in the YAML), then:\n'
          f'  snobal-box run {os.path.join(args.dir, "my_site.yaml")}')


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog='snobal-box',
        description='Run the snobal point snow model from a YAML config.')
    sub = ap.add_subparsers(dest='cmd', required=True)

    r = sub.add_parser('run', help='run baseline and scenarios, write output')
    r.add_argument('config', nargs='+', help='YAML config file(s)')
    r.add_argument('-s', '--scenario', action='append',
                   help='run only this scenario label (repeatable)')
    r.add_argument('--combined', help='also write all summaries to this CSV')
    r.add_argument('-q', '--quiet', action='store_true')
    r.set_defaults(fn=_run)

    c = sub.add_parser('check', help='read and QC the forcing; no model run')
    c.add_argument('config', nargs='+')
    c.set_defaults(fn=_check)

    i = sub.add_parser('init', help='write a template config and CSV')
    i.add_argument('dir', nargs='?', default='my_site')
    i.add_argument('--force', action='store_true')
    i.set_defaults(fn=_init)

    args = ap.parse_args(argv)
    warnings.simplefilter('default')
    try:
        args.fn(args)
    except (ValueError, FileNotFoundError) as e:
        print(f'error: {e}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
