.. :changelog:

History
=======

0.3.1 (unreleased)
------------------

* Fix issues #108 #109 and merges #110 #120
* Remove requirements.txt, add docs/rtd-pip-requirements, update .travis.yaml
* Add support for running at sites where data is locally visible (e.g., nersc)
* Add support for reading spArc and spFlat calibration files.
* Add CalibrationTutorial notebook.
* Add bosssky script.
* Fix some issues with python3 and astropy < 2.

0.3.0 (2017-12-17)
------------------

* Update for python 3.x

0.2.8 (2015-10-19)
------------------

* Fix issues #29 #103 and addresses #106
* Add support for reading raw image data and plug map files.
* Add warnings about not using numpy 1.10.0 or 1.10.1.

0.2.7 (2015-09-28)
------------------

* Fix issues #92 #94 #96 #97 #100
* Add support for reading per-exposure flux calibration and correction vectors.
* Add plot functions for per-fiber data vs fiber number or focal-plane position.
* Add a plug_map attribute to spPlate, spFrame, spCFrame.
* FrameFile infers the spectrograph index and whether flux calibration has been applied.
* bossdata infers MJD when possible.
* bossplot option "--camera" renamed to "--band".

0.2.6 (2015-08-05)
------------------

* Fix issues #67 #74 #86
* The ``camera`` arg to ``SpecFile.get_valid_data`` (and related methods) should now be ``b1``, ``b2``, ``r1``, ``r2`` instead of ``blue`` or ``red``.
* New options for the ``get_valid_data`` methods: ``use_ivar``, ``use_loglam``, ``fiducial_grid``.

0.2.5 (2015-07-06)
------------------

* Fix issues #27 #28 #63 #64 #68
* New command-line options include:

 * bossplot: --platefile, --flux-range, --wlen-range, --wdisp-range, --label-pos, --no-grid, --show-invalid
 * bossfetch: --platefile

* Adds support for spPlate files and platelist metadata.
* Adds command-line options to customize bossplot axes, add labels and grids, and display invalid data.
* General documentation cleanup.
* Better error handling in bossplot.

0.2.4 (2015-06-29)
------------------

* Fix issues #11 #36 #41 #43 #45 #50
* New command-line options include:

 * bossfetch: --plate-name, --mjd-name, --fiber-name
 * bosscatalog: --quasar-catalog, --quasar-catalog-name

* The main new functionality is support for querying the quasar catalog, using different data sources, and built-in defaults for any of the four environment variables that is not set.

0.2.3 (2015-06-22)
------------------

* Fix issues #2 #10 #16 #18 #19 #21 #24
* New command-line options include:

 * bossfetch: --globus, --dry-run
 * bossplot: --save-data
 * bossquery: --sort

* The main new library functionality is support for using wavelengths and dispersions encoded as "trace sets" in spFrame files via :class:`bossdata.plate.TraceSet`.

0.2.2 (2015-06-15)
------------------

* Really fix issues #9 #13.
* Add support for finding and fetching spFrame and spCFrame files (#17).

0.2.1 (2015-06-13)
------------------

* Fix issues #9 #12 #13

0.2.0 (2015-06-09)
------------------

* Fix issues #3 #5 #6
* Add support for accessing subtracted sky flux to the `spec` module and `bossplot` script.
* This version breaks backwards compatiblity with 0.1.0 since the previous `$BOSS_SAS_ROOT` environment variable is now named `$BOSS_SAS_PATH` and has the instrument name (usually `boss`) appended.
* bash users can update by replacing `export BOSS_SAS_ROOT=/sas/dr12` with `export BOSS_SAS_PATH=/sas/dr12/boss` in their `.bashrc` file.

0.1.0 (2015-05-24)
------------------

* First release on PyPI.
