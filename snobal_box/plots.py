"""Figures for a RunResult. Needs matplotlib (installed with the package)."""

import os

import matplotlib
import matplotlib.dates as mdates

__all__ = ['plot_swe_depth', 'plot_water', 'plot_energy', 'save_figures']


def _plt():
    import matplotlib.pyplot as plt
    return plt


def _tidy(ax, ylabel, title=None):
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title, loc='left', fontsize=11)
    ax.grid(alpha=0.25, lw=0.5)
    ax.xaxis.set_major_locator(mdates.MonthLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)


def plot_swe_depth(outputs, title='', observed=None, figsize=(11, 6.5)):
    """``outputs`` is one output frame or {label: frame}. ``observed`` is a
    dict that may hold 'depth' (m) and 'swe' (mm) Series."""
    plt = _plt()
    runs = outputs if isinstance(outputs, dict) else {'model': outputs}
    observed = observed or {}
    fig, axes = plt.subplots(2, 1, figsize=figsize, sharex=True)
    for label, o in runs.items():
        lw = 1.8 if label == 'baseline' else 1.2
        axes[0].plot(o.index, o['m_s'], lw=lw, label=label)
        axes[1].plot(o.index, o['z_s'], lw=lw, label=label)
    if 'swe' in observed:
        s = observed['swe']
        axes[0].plot(s.index, s.values, 'k.', ms=3, label='observed')
    if 'depth' in observed:
        s = observed['depth']
        axes[1].plot(s.index, s.values, 'k.', ms=3, label='observed')
    _tidy(axes[0], 'SWE (mm)', f'{title} snow water equivalent'.strip())
    _tidy(axes[1], 'depth (m)', 'Snow depth')
    if len(runs) > 1 or observed:
        axes[0].legend(fontsize=8, frameon=False)
    fig.tight_layout()
    return fig


def plot_water(out, title='', figsize=(11, 6.5)):
    plt = _plt()
    fig, axes = plt.subplots(2, 1, figsize=figsize, sharex=True)
    axes[0].plot(out.index, out['m_s'], color='#2b6ca3', lw=1.5)
    axes[0].fill_between(out.index, 0, out['m_s'], color='#2b6ca3',
                         alpha=0.15)
    _tidy(axes[0], 'SWE (mm)', f'{title} snowpack'.strip())
    axes[1].plot(out.index, out['ro_pred_sum'].cumsum(), lw=1.6,
                 color='#2b6ca3', label='snow water input')
    axes[1].plot(out.index, out['melt_sum'].cumsum(), lw=1.3,
                 color='#F9A620', label='melt (includes refrozen)')
    axes[1].plot(out.index, -out['E_s_sum'].cumsum(), lw=1.3,
                 color='#8c6d4f', label='sublimation/evaporation')
    _tidy(axes[1], 'cumulative (mm)', 'Water leaving the snowpack')
    axes[1].legend(fontsize=9, frameon=False)
    fig.tight_layout()
    return fig


def plot_energy(out, title='', window=24 * 7, figsize=(11, 5)):
    plt = _plt()
    terms = [('R_n_bar', 'net radiation', '#d62728'),
             ('H_bar', 'sensible', '#ff7f0e'),
             ('L_v_E_bar', 'latent', '#1f77b4'),
             ('G_bar', 'ground', '#8c564b'),
             ('delta_Q_bar', 'net', '#000000')]
    has_snow = out['m_s'] > 0
    fig, ax = plt.subplots(figsize=figsize)
    for col, label, color in terms:
        sm = out[col].where(has_snow).rolling(
            window, center=True, min_periods=max(1, window // 4)).mean()
        ax.plot(out.index, sm.where(has_snow), color=color,
                lw=1.8 if col == 'delta_Q_bar' else 1.2, label=label)
    ax.axhline(0, color='k', lw=0.6, ls=':')
    _tidy(ax, 'W m$^{-2}$', f'{title} energy balance, {window // 24}-day '
          'mean, positive into the pack'.strip())
    ax.legend(fontsize=9, frameon=False, ncol=5)
    fig.tight_layout()
    return fig


def save_figures(result, outdir, dpi=130):
    matplotlib.use('Agg', force=False)
    plt = _plt()
    name = result.config['name']
    base = result.outputs.get('baseline',
                              next(iter(result.outputs.values())))
    figs = {
        'fig_swe_depth': plot_swe_depth(result.outputs, name,
                                        result.observations),
        'fig_water': plot_water(base, name),
        'fig_energy': plot_energy(base, name),
    }
    paths = {}
    for key, fig in figs.items():
        p = os.path.join(outdir, f'{key[4:]}.png')
        fig.savefig(p, dpi=dpi, bbox_inches='tight')
        plt.close(fig)
        paths[key] = p
    return paths
