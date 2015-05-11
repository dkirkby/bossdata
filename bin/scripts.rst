Executable scripts
==================

Before running any scripts, you will normally need to establish your local configuration with some environment variables:

 * BOSS_DATA_ROOT: The top-level directory where all downloaded data files will be locally mirrored. Make sure there is enough space here for the files you plan to use locally. You might want to exclude this directory from your backups since it can get large and can be easily recreated.
 * BOSS_DATA_URL: The top-level URL for downloading BOSS data. This should normally be `http://dr12.sdss3.org` for publicly accessible data.
 * BOSS_SAS_ROOT: This identifies the data release that you want to work with. This should normally be /sas/dr12, which corresponds to the final BOSS Data Release 12.
 * BOSS_REDUX_VERSION: This is the pipeline reconstruction version that you want to work with. This should normally be either v5_7_0, which corresponds to the final processing of all BOSS data, or v5_7_2 for SEQUELS data.

As a sanity check of your configuration, `$BOSS_DATA_ROOT` should be an existing directory where you have write access (and plenty of free space), and `$BOSS_DATA_URL/$BOSS_SAS_ROOT/boss/spectro/redux/$BOSS_REDUX_VERSION/` should be a valid URL that displays a directory listing in any browser.

For complete documentation on the command-line options of any script use the `--help` option, for example::

    bossquery --help

bossquery
---------

Query the meta data for BOSS observations. For example::

    bossquery --what PLATE,MJD,FIBER,PLUG_RA,PLUG_DEC,Z --where 'OBJTYPE="QSO"' --save qso.dat

The save option supports `many different output formats <http://astropy.readthedocs.org/en/latest/io/unified.html#built-in-table-readers-writers>`_ that are automatically selected based on the file extension.  In addition, this program automatically maps the `.dat` and `.txt` extensions to the `ascii` format.

This command uses an sqlite3 database of metadata that will be created if necessary. By default, the "lite" version database will be used, which provides faster queries and a smaller database file.  However, the full `spAll data model <http://dr12.sdss3.org/datamodel/files/BOSS_SPECTRO_REDUX/RUN2D/spAll.html>`_ is also available with the `--full` option (resulting in slower queries and a larger database file).  The "lite" and "full" databases are separate files based on different downloads. Once either has been created the first time, it will be immediately available for future queries.

bossfetch
---------

Fetch BOSS data files containing the spectra of specified observations and mirror them locally. For example::

    bossfetch --verbose qso.dat

By default, the "lite" format of each spectrum data file is downloaded, which is sufficient for many purposes and signficantly (about 8x) smaller. The "lite" format contains HDUs 0-3 of the `full spectrum data file <http://dr12.sdss3.org/datamodel/files/BOSS_SPECTRO_REDUX/RUN2D/spectra/PLATE4/spec.html>`_ and does not include the spectra of individual exposures.  To download the full files instead, use the `--full` option. Both types of files can co-exist in your local mirror.

The `--verbose` option displays a progress bar showing the fraction of files already locally available. Any files that were previously fetched will not be downloaded again so it is safe and efficient to run `bossfetch` for overlapping lists of observations.  Note that the progress bar may appear to update unevenly if some files are already mirrored and others need to be downloaded.

Each data file download is streamed to a temporary files with `.downloading` appended to their name then renamed to remove this extension after the download completes normally. If a download is interrupted or fails for some reason, the partially downloaded file will remain in the local mirror.  Re-running a `bossfetch` command will automatically re-download any partially downloaded file.
