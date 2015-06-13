.. :changelog:

History
-------

0.1.0 (2015-05-24)
------------------

* First release on PyPI.

0.2.0 (2015-06-09)
------------------

* Fix issues #3 #5 #6.
* Add support for accessing subtracted sky flux to the `spec` module and `bossplot` script.
* This version breaks backwards compatiblity with 0.1.0 since the previous `$BOSS_SAS_ROOT` environment variable is now named `$BOSS_SAS_PATH` and has the instrument name (usually `boss`) appended.
* bash users can update by replacing `export BOSS_SAS_ROOT=/sas/dr12` with `export BOSS_SAS_PATH=/sas/dr12/boss` in their `.bashrc` file.

0.2.1 (2015-06-13)
------------------

* Fix issues #9 #12 #13.
