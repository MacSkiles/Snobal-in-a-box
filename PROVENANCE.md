# Provenance

`snobal_box/_pysnobal/` is the physics core of **PySnobal 0.2.3**
(scotthavens/pysnobal, released 2020-07-29), distributed under **CC0 1.0**
as a work of the United States government. PySnobal is a Python port of the
snobal point model (Marks et al. 1999, Hydrological Processes 13:1935-1959).

Retained: `point/` (the physics) and `core/constants.py`.
Removed: the config/CLI layer, netCDF output and the gridded `spatial`
module, which drops the inicheck, pytz and setuptools_scm dependencies and
the `numpy < 1.19` pin.

Local modifications to the vendored files: import paths rewritten to
`snobal_box._pysnobal`. No physics changes.

Run-time substitution (not an edit of the vendored code): `snobal_box.runner`
wraps `libsnobal.hle1` so that, when the turbulent-transfer iteration fails to
converge, neutral-stability fluxes are used for that step instead of aborting
the run. Every substitution is counted in the run summary
(`neutral_stability_steps`).

Everything outside `_pysnobal/` is original to this package.

TODO before publishing: confirm the "no physics changes" line against a
diff with upstream 0.2.3.
