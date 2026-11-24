"""Vendored physics core of PySnobal 0.2.3.

Upstream: scotthavens/pysnobal, released 2020-07-29 under CC0 1.0
(a work of the United States government, in the public domain within the US,
with copyright and related rights waived worldwide).

Only `point` (the physics) and `core.constants` are retained. The original
config, CLI and netCDF layer and the gridded `spatial` module are removed,
which also removes the inicheck, pytz and setuptools_scm dependencies and the
`numpy < 1.19` pin that makes upstream uninstallable on modern Python.

Local modifications are recorded in PROVENANCE.md at the repository root.
Nothing in this subpackage is part of the public API; import from
`snobal_box` instead.
"""
