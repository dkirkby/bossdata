Executable scripts
==================

Before running any scripts, you will normally need to establish your local configuration with some environment variables:

* BOSS_LOCAL_ROOT: The top-level directory where all downloaded data files will be locally mirrored. Make sure there is enough space here for the files you plan to use locally. You might want to exclude this directory from your backups since it can get large and can be easily recreated.
* BOSS_DATA_URL: The top-level URL for downloading BOSS data. This should normally be `http://dr12.sdss3.org` for publicly accessible data.
* BOSS_SAS_PATH: This identifies the top-level path under SAS of the data you want to work with, and should normally begin with "/sas". For the final BOSS Data Release 12, use "/sas/dr12/boss".
* BOSS_REDUX_VERSION: This is the pipeline reconstruction version that you want to work with. This should normally be either v5_7_0, which corresponds to the final processing of all BOSS data, or v5_7_2 for SEQUELS data.

As a sanity check of your configuration, the following shell command::

    echo $BOSS_DATA_URL/$BOSS_SAS_PATH/boss/spectro/redux/$BOSS_REDUX_VERSION/

should print a valid URL that displays a directory listing in any browser, without requiring any authentication. For bash users, the following lines added to your `~/.bashrc` file should be a good starting point::

    export BOSS_LOCAL_ROOT=...some suitable local path...
    export BOSS_DATA_URL=http://dr12.sdss3.org
    export BOSS_SAS_PATH=/sas/dr12/boss
    export BOSS_REDUX_VERSION=v5_7_0

For complete documentation on the command-line options of any script use the `--help` option, for example::

    bossquery --help

.. _bossquery:

bossquery
---------

Query the meta data for BOSS observations. For example::

    bossquery --what PLATE,MJD,FIBER,PLUG_RA,PLUG_DEC,Z --where 'OBJTYPE="QSO"' --save qso.dat

The save option supports `many different output formats <http://astropy.readthedocs.org/en/latest/io/unified.html#built-in-table-readers-writers>`_ that are automatically selected based on the file extension.  In addition, this program automatically maps the `.dat` and `.txt` extensions to the `ascii` format.

This command uses an sqlite3 database of metadata that will be created if necessary. By default, the "lite" version database will be used, which provides faster queries and a smaller database file.  However, the full `spAll data model <http://dr12.sdss3.org/datamodel/files/BOSS_SPECTRO_REDUX/RUN2D/spAll.html>`_ is also available with the `--full` option (resulting in slower queries and a larger database file).  The "lite" and "full" databases are separate files based on different downloads. Once either has been created the first time, it will be immediately available for future queries.  Note that it can take a while to create the initial database file: allow about 30 minutes for either version. Once the database has been created, you can safely delete the downloaded source file if you are short on disk space.

The columns in the lite database are a subset of those in the full database but the values are not numerically identical between them because they are truncated in the text file used to generate the lite database. However, the level of these truncation errors should be insigificant for any science applications.

There are some minor inconsistencies between the data models of the lite and full versions of the meta data provided by BOSS.  In particular, the lite format uses the name `FIBER` while the full version uses `FIBERID`. We resolve this by consistently using the shorter form `FIBER` in both SQL databases.  Also, the full format includes columns that are themselves arrays. One of these, `MODELFUX(5)`, is included in the lite format using names `MODELFLUX0...MODELFUX4`. We normalize the mapping of array columns to scalar SQL columns using the syntax `COLNAME_I` for element [i] of a 1D array and `COLNAME_I_J` for element [i,j] of a 2D array, with indices starting from zero. This means, for example, that `MODELFLUX(5)` values are consistently named `MODELFLUX_0...MODELFLUX_4` in both SQL databases.

.. _bossfetch:

bossfetch
---------

Fetch BOSS data files containing the spectra of specified observations and mirror them locally. For example::

    bossfetch --verbose qso.dat

Fetched files will be placed under `$BOSS_LOCAL_ROOT` with paths that exactly match the URLs they are downloaded from with the prefix substitution::

    $BOSS_DATA_URL => $BOSS_LOCAL_ROOT

For example, with the default configuration given above, the file at::

    http://dr12.sdss3.org/sas/dr12/boss/spectro/redux/v5_7_0/spectra/lite/3586/spec-3586-55181-0190.fits

would be downloaded to::

    $BOSS_LOCAL_ROOT/sas/dr12/boss/spectro/redux/v5_7_0/spectra/lite/3586/spec-3586-55181-0190.fits

By default, the "lite" format of each spectrum data file is downloaded, which is sufficient for many purposes and signficantly (about 8x) smaller. The "lite" format contains HDUs 0-3 of the `full spectrum data file <http://dr12.sdss3.org/datamodel/files/BOSS_SPECTRO_REDUX/RUN2D/spectra/PLATE4/spec.html>`_ and does not include the spectra of individual exposures.  To download the full files instead, use the `--full` option. Both types of files can co-exist in your local mirror.

The `--verbose` option displays a progress bar showing the fraction of files already locally available. Any files that were previously fetched will not be downloaded again so it is safe and efficient to run `bossfetch` for overlapping lists of observations.  Note that the progress bar may appear to update unevenly if some files are already mirrored and others need to be downloaded.

Each data file download is streamed to a temporary files with `.downloading` appended to their name then renamed to remove this extension after the download completes normally. If a download is interrupted or fails for some reason, the partially downloaded file will remain in the local mirror.  Re-running a `bossfetch` command will automatically re-download any partially downloaded file.

By default, downloading is split between two parallel subprocesses but you can change this with the
`--nproc` option.  For downloading "lite" files, using more than 2 subprocesses will probably not
improve the overall performance.

.. _bossplot:

bossplot
--------

Plot the spectrum of a single BOSS observation, identified by its PLATE, MJD of the observation, and the FIBER that was assigned to the target whose spectrum you want to plot. For example (these are the defaults if you omit any parameters)::

    bosplot --plate 6641 --mjd 56383 --fiber 30

This should open a new window containing the plot that you will need to close in order to exit the program.  To also save your plot, add the `--save` option with a filename that has a standard graphics format extension (pdf,png,...).  To save plots directly without displaying them, also use the `--no-display` option.

The `bossplot` command will automatically download the appropriate data file if necessary.

This script uses the `matplotlib` python library, which is not required for the `bossdata` package and therefore not automatically installed.
