.. :changelog:

History
=======

0.1.0 (2015-05-24)
------------------

* First release on PyPI.

0.2.0 (2015-06-09)
------------------

* Fix issues #3 #5 #6
* Add support for accessing subtracted sky flux to the `spec` module and `bossplot` script.
* This version breaks backwards compatiblity with 0.1.0 since the previous `$BOSS_SAS_ROOT` environment variable is now named `$BOSS_SAS_PATH` and has the instrument name (usually `boss`) appended.
* bash users can update by replacing `export BOSS_SAS_ROOT=/sas/dr12` with `export BOSS_SAS_PATH=/sas/dr12/boss` in their `.bashrc` file.

0.2.1 (2015-06-13)
------------------

* Fix issues #9 #12 #13

0.2.2 (2015-06-15)
------------------

* Really fix issues #9 #13.
* Add support for finding and fetching spFrame and spCFrame files (#17).

0.2.3 (2015-06-22)
------------------

* Fix issues #2 #10 #16 #18 #19 #21 #24
* New command-line options include:

 * bossfetch: --globus, --dry-run
 * bossplot: --save-data
 * bossquery: --sort

* The main new library functionality is support for using wavelengths and dispersions encoded as "trace sets" in spFrame files via :class:`bossdata.plate.TraceSet`.

0.2.4 (2015-06-29)
------------------

* Fix issues #11 #36 #41 #43 #45 #50
* New command-line options include:

 * bossfetch: --plate-name, --mjd-name, --fiber-name
 * bosscatalog: --quasar-catalog, --quasar-catalog-name

* The main new functionality is support for querying the quasar catalog, using different data sources, and built-in defaults for any of the four environment variables that is not set.

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
