"""One Site holds everything needed to run one location."""

import os
from dataclasses import dataclass
from typing import Optional

from snobal_box.io import (read_forcing, read_initial_snow,
                           read_measurement_heights, read_reference_output)
from snobal_box.phase import repair_precip  # noqa: F401
from snobal_box.perturb import perturb
from snobal_box.model import run

__all__ = ['Site']

@dataclass
class Site:
    """Everything needed to run one site. Point `root` at wherever the files
    live -- a local folder, a mounted Google Drive path, anything.

    Filenames are resolved relative to `root` unless they are absolute.
    """
    name: str
    root: str
    data_file: str
    precip_file: str
    snow_file: str
    height_file: str
    elevation: float
    model_file: Optional[str] = None
    start_date: str = '2000-10-01'
    tstep_minutes: int = 60
    temps_in_C: bool = True
    max_h2o_vol: float = 0.01
    max_z_s_0: float = 0.25
    start_record: int = 0      # skip this many leading forcing records
    n_records: Optional[int] = None   # truncate after this many
    notes: str = ''

    def path(self, filename):
        if filename is None:
            return None
        return filename if os.path.isabs(filename) \
            else os.path.join(self.root, filename)

    def check(self):
        """Verify every referenced file exists before doing any work."""
        missing = []
        for attr in ['data_file', 'precip_file', 'snow_file', 'height_file']:
            p = self.path(getattr(self, attr))
            if not os.path.exists(p):
                missing.append(f'  {attr}: {p}')
        if missing:
            raise FileNotFoundError(
                f'{self.name}: missing input file(s)\n' + '\n'.join(missing) +
                f'\n\nFiles actually present in {self.root!r}:\n  ' +
                '\n  '.join(sorted(os.listdir(self.root))
                            if os.path.isdir(self.root)
                            else ['<root is not a directory>']))
        return True

    def load(self, repair=True, verbose=True):
        """Returns (forcing, snow_prop, meas_heights, reference_or_None)."""
        self.check()
        mh = read_measurement_heights(self.path(self.height_file))
        snow = read_initial_snow(self.path(self.snow_file), z_0=mh['z_0'],
                                 temps_in_C=self.temps_in_C)
        forcing = read_forcing(self.path(self.data_file),
                               self.path(self.precip_file),
                               start_date=self.start_date,
                               tstep_minutes=self.tstep_minutes,
                               temps_in_C=self.temps_in_C)
        if self.start_record or self.n_records:
            end = None if self.n_records is None \
                else self.start_record + self.n_records
            forcing = forcing.iloc[self.start_record:end]
        if repair:
            forcing = repair_precip(forcing, verbose=verbose)
        ref = read_reference_output(self.path(self.model_file))
        if verbose:
            print(f'{self.name}: {len(forcing)} records, '
                  f'{forcing.index[0].date()} to {forcing.index[-1].date()} '
                  f'({len(forcing) / 24:.0f} days)')
            print(f'  initial SWE {snow["z_s"] * snow["rho"]:.1f} mm '
                  f'(depth {snow["z_s"]:.2f} m, density {snow["rho"]:.0f})')
            print(f'  total precip {forcing["m_pp"].sum():.1f} mm in '
                  f'{(forcing["m_pp"] > 0).sum()} events')
            print(f'  reference output: '
                  f'{"present" if ref is not None else "none"}')
        return forcing, snow, mh, ref

    def run(self, forcing=None, tstep_info=None, **perturbations):
        """Load (if needed), optionally perturb, and run."""
        if forcing is None:
            forcing, snow, mh, _ = self.load(verbose=False)
        else:
            _, snow, mh, _ = self.load(verbose=False)
        if perturbations:
            forcing = perturb(forcing, **perturbations)
        return run(forcing, snow, mh, elevation=self.elevation,
                   max_h2o_vol=self.max_h2o_vol, max_z_s_0=self.max_z_s_0,
                   tstep_info=tstep_info)
