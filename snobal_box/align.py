"""Find or check the start time of an untimestamped forcing record.

IPW ``.input`` files carry an hourly counter, not times, so their start time
is an assertion. If a timestamped record of the same place exists (an HRRR
extraction, a nearby station), the diurnal cycle of air temperature pins the
start to the hour: slide one series along the other and keep the lag with the
highest correlation.
"""

import numpy as np
import pandas as pd

__all__ = ['find_start']


def find_start(series, reference, guess=None, search_days=45,
               min_overlap_hours=24 * 30):
    """Best start time for ``series`` (an array of hourly values, no times)
    against ``reference`` (an hourly Series with a tz-aware index).

    Candidate starts are every hour within ``search_days`` of ``guess``
    (default: the reference's first time). Returns ``(start, r, overlap)``
    with ``start`` in the reference's timezone.

    Correlation is computed on anomalies from a 24 h rolling mean, so the
    match comes from timing of the diurnal cycle rather than the seasonal
    trend.
    """
    ref = reference.asfreq('h')
    a = np.asarray(series, dtype=float)
    a = a - pd.Series(a).rolling(24, center=True, min_periods=1).mean()\
        .to_numpy()
    b = ref - ref.rolling(24, center=True, min_periods=1).mean()
    b = b.to_numpy()

    guess = ref.index[0] if guess is None else pd.Timestamp(guess)
    if guess.tzinfo is None:
        guess = guess.tz_localize(ref.index.tz)
    i0 = int(round((guess - ref.index[0]) / pd.Timedelta('1h')))
    span = int(search_days * 24)

    best = (None, -np.inf, 0)
    for k in range(i0 - span, i0 + span + 1):
        # series[j] aligns with ref[k + j]
        j0 = max(0, -k)
        j1 = min(len(a), len(b) - k)
        if j1 - j0 < min_overlap_hours:
            continue
        x, y = a[j0:j1], b[k + j0:k + j1]
        ok = np.isfinite(x) & np.isfinite(y)
        if ok.sum() < min_overlap_hours:
            continue
        r = np.corrcoef(x[ok], y[ok])[0, 1]
        if r > best[1]:
            best = (k, r, int(ok.sum()))
    if best[0] is None:
        raise ValueError('no lag with enough overlap; widen search_days or '
                         'check that the reference covers the record')
    k, r, n = best
    start = ref.index[0] + pd.Timedelta(hours=k)
    return start, float(r), n
