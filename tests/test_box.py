import os
import warnings

import numpy as np
import pandas as pd
import pytest

import snobal_box as sb
from snobal_box.forcing import read_ipw_forcing, read_csv_forcing
from snobal_box.precip import apply_precip_scheme, PHASE_SCHEMES
from snobal_box.runner import scenario_forcing

HERE = os.path.dirname(__file__)
ROOT = os.path.dirname(HERE)
TEMPLATE_CSV = os.path.join(ROOT, 'snobal_box', 'templates',
                            'forcing_template.csv')


def _write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return str(p)


def test_ipw_empty_field_is_repaired_not_shifted(tmp_path):
    data = _write(tmp_path, 'd.input', '0.1 288 7.5 513 0.5 -2.5\n' * 3)
    precip = _write(tmp_path, 'p.input',
                    '0  0.0 0.0 5.24\n1 1.0 1.0 100.0 -3.0\n2 0.0 0.0 0.0 4\n')
    start = pd.Timestamp('2023-10-01', tz='Etc/GMT+7')
    with pytest.warns(UserWarning, match='empty field'):
        f = read_ipw_forcing(data, precip, start)
    assert f['m_pp'].tolist() == [0.0, 1.0, 0.0]
    assert f['rho_snow'].tolist() == [0.0, 100.0, 0.0]   # not a temperature
    assert f.index[0] == start


def test_ipw_bad_row_is_an_error(tmp_path):
    data = _write(tmp_path, 'd.input', '0.1 288 7.5 513\n')
    precip = _write(tmp_path, 'p.input', '0 0 0 0 0\n')
    with pytest.raises(ValueError, match='line 1'):
        read_ipw_forcing(data, precip, pd.Timestamp('2023-10-01', tz='UTC'))


def test_csv_template_reads_and_gaps(tmp_path):
    f = read_csv_forcing(TEMPLATE_CSV)
    assert len(f) == 48 and str(f.index.tz) == 'UTC'
    d = pd.read_csv(TEMPLATE_CSV).drop(index=[5, 6])
    p = tmp_path / 'gappy.csv'
    d.to_csv(p, index=False)
    with pytest.warns(UserWarning, match='filled 2'):
        g = read_csv_forcing(str(p))
    assert len(g) == 48
    d = pd.read_csv(TEMPLATE_CSV).drop(index=range(5, 15))
    d.to_csv(p, index=False)
    with pytest.raises(ValueError, match='longest run'):
        read_csv_forcing(str(p))


def _frame(T_C):
    T = np.asarray(T_C, dtype=float) + 273.16
    n = len(T)
    return pd.DataFrame({'S_n': 0.0, 'I_lw': 250.0, 'T_a': T, 'e_a': 400.0,
                         'u': 1.0, 'T_g': 270.66, 'm_pp': 1.0, 'T_pp': T},
                        index=pd.date_range('2024-01-01', periods=n,
                                            freq='h', tz='UTC'))


@pytest.mark.parametrize('scheme', sorted(PHASE_SCHEMES))
def test_phase_schemes_are_monotone_and_bounded(scheme):
    f = apply_precip_scheme(_frame(np.linspace(-10, 10, 81)),
                            {'phase': scheme, 'density': 'hedstrom_pomeroy'})
    p = f['percent_snow'].to_numpy()
    assert (p >= 0).all() and (p <= 1).all()
    assert (np.diff(p) <= 1e-12).all()        # colder never less snow
    assert p[0] > 0.95 and p[-1] < 0.05
    assert ((f['rho_snow'] > 0) == (f['percent_snow'] > 0)).all()


def test_threshold_shift_positive_means_more_snow():
    base = apply_precip_scheme(_frame([0.5]), {'phase': 'threshold'})
    more = apply_precip_scheme(_frame([0.5]),
                               {'phase': 'threshold', 'threshold_shift': 1})
    assert base['percent_snow'].iloc[0] == 0 and \
        more['percent_snow'].iloc[0] == 1


@pytest.mark.parametrize('with_tpp', [True, False])
def test_energy_only_scenario_keeps_phase(with_tpp):
    raw = _frame([-0.8] * 5)
    if not with_tpp:
        raw = raw.drop(columns='T_pp')
    pc = {'phase': 'table', 'density': 'table'}
    warm, _ = scenario_forcing(raw, pc, T_offset=2.0)
    held, _ = scenario_forcing(raw, pc, T_offset=2.0, rebuild_phase=False)
    base, _ = scenario_forcing(raw, pc)
    assert warm['percent_snow'].iloc[0] < base['percent_snow'].iloc[0]
    assert held['percent_snow'].iloc[0] == base['percent_snow'].iloc[0]
    assert (held['T_a'] == warm['T_a']).all()
    assert np.allclose(held['T_pp'], warm['T_pp'])


def test_find_start_recovers_known_offset():
    t = pd.date_range('2021-10-01', periods=24 * 90, freq='h', tz='UTC')
    rng = np.random.default_rng(0)
    ref = pd.Series(5 * np.sin(2 * np.pi * (t.hour - 15) / 24)
                    + rng.normal(0, 1, len(t)), index=t)
    series = ref.to_numpy()[31:31 + 24 * 60] + rng.normal(0, 0.5, 24 * 60)
    start, r, _ = sb.find_start(series, ref, search_days=5)
    assert start == t[31] and r > 0.8


def test_config_errors_are_listed():
    with pytest.raises(ValueError) as e:
        sb.load_config({'forcing': {'format': 'ipw', 'data_file': 'a',
                                    'precip_file': 'b'}})
    msg = str(e.value)
    for needle in ('name', 'elevation', 'heights', 'forcing.start'):
        assert needle in msg


def test_short_csv_run_closes_mass(tmp_path):
    cfg = {
        'name': 'template', 'site': {'elevation_m': 2670,
                                     'heights': {'z_u': 5, 'z_t': 5,
                                                 'z_0': 0.003, 'z_g': 0.15},
                                     'initial_snow': {'depth_m': 0.5,
                                                      'density': 300}},
        'forcing': {'format': 'csv', 'file': TEMPLATE_CSV},
        'precip': {'phase': 'linear', 'density': 'hedstrom_pomeroy'},
        'scenarios': {'+2C': {'T_offset': 2}},
        'output': {'dir': str(tmp_path), 'figures': False},
    }
    res = sb.run_config(cfg, verbose=False)
    s = res.summary
    assert list(s.index) == ['baseline', '+2C']
    assert (s['closure_residual_mm'].abs() < 0.5).all()
    assert os.path.exists(res.paths['summary'])
