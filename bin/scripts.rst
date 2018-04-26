Executable scripts
==================

For complete documentation on the command-line options of any script use the `--help` option, for example::

    bossquery --help

You will normally want to configure ``bossdata`` by setting some :doc:`environment variables </envvar>`.

.. _bossquery:

bossquery
---------

Query the meta data for BOSS observations. For example (the initial setup might take a few minutes if you have never used this command before)::

    bossquery --what PLATE,MJD,FIBER,PLUG_RA,PLUG_DEC,Z --where 'OBJTYPE="QSO"' --sort Z --save qso.dat

The `--save` option supports `many different output formats <http://astropy.readthedocs.org/en/latest/io/unified.html#built-in-table-readers-writers>`_ that are automatically selected based on the file extension.  In addition, this program automatically maps the `.dat` and `.txt` extensions to the `ascii` format.

The `--what`, `--where` and `--sort` options all use SQL syntax (these are in fact substituted into a SQL string).

* `--what` takes a comma separated list of column names (like SQL SELECT) and defaults to PLATE,MJD,FIBER::

    --what PLATE,MJD,FIBER,PLUG_RA,PLUG_DEC,Z

* `--where` takes a SQL 'WHERE' string::

    --where '(OBJTYPE="QSO" and Z > 0.1) or CLASS="QSO"'

* `--sort` takes a list of columns with optional DESC keyword following columns to reverse their order (a la SQL ORDER BY)::

    --sort 'CLASS, Z DESC'

This command uses an sqlite3 database of metadata that will be created if necessary. By default, the "lite" version database will be used, which provides faster queries and a smaller database file.  However, the full `spAll data model <http://dr12.sdss3.org/datamodel/files/BOSS_SPECTRO_REDUX/RUN2D/spAll.html>`_ is also available with the `--full` option (resulting in slower queries and a larger database file).  The "lite" and "full" databases are separate files based on different downloads. Once either has been created the first time, it will be immediately available for future queries.  Note that it can take a while to create the initial database file: allow about 30 minutes for either version. Once the database has been created, you can safely delete the downloaded source file if you are short on disk space.

The columns in the lite database are a subset of those in the full database but the values are not numerically identical between them because they are truncated in the text file used to generate the lite database. However, the level of these truncation errors should be insignificant for any science applications.

There are some minor inconsistencies between the data models of the lite and full versions of the meta data provided by BOSS.  In particular, the lite format uses the name `FIBER` while the full version uses `FIBERID`. We resolve this by consistently using the shorter form `FIBER` in both SQL databases.  Also, the full format includes columns that are themselves arrays. One of these, `MODELFUX(5)`, is included in the lite format using names `MODELFLUX0...MODELFUX4`. We normalize the mapping of array columns to scalar SQL columns using the syntax `COLNAME_I` for element [i] of a 1D array and `COLNAME_I_J` for element [i,j] of a 2D array, with indices starting from zero. This means, for example, that `MODELFLUX(5)` values are consistently named `MODELFLUX_0...MODELFLUX_4` in both SQL databases.

In the case where a query is made without specifying `--full` but the lite database file is not present, an attempt will be made to use the full database.  If neither DB files are present the same logic is applied to the catalog files.  If present, the lite catalog file will be parsed and the lite DB created; if that is not present, the full catalog file will be parsed and the full DB created.  Only after exhausting these options will a download (of the lite DB) file be attempted.

Note that specifying `--full` will only (and always) use the full DB or catalog file.

The `--quasar-catalog` option can be used to query the `BOSS quasar catalog <http://www.sdss.org/dr12/algorithms/boss-dr12-quasar-catalog/>`_ instead of spAll. By default, the current version of the catalog will be used; use the `--quasar-catalog-name` option to specify an earlier version.

The ```--platelist`` option can be used to query the `BOSS plate list database <http://data.sdss3.org/datamodel/files/BOSS_SPECTRO_REDUX/RUN2D/platelist.html>`_ instead of spAll.

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

By default, the "lite" format of each spectrum data file is downloaded, which is sufficient for many purposes and signficantly (about 8x) smaller. The "lite" format contains HDUs 0-3 of the `full spectrum data file <http://dr12.sdss3.org/datamodel/files/BOSS_SPECTRO_REDUX/RUN2D/spectra/PLATE4/spec.html>`_ and does not include the spectra of individual exposures.  To download the full files instead, use the ``--full`` option. Both types of files can co-exist in your local mirror. You can also load the plate ``spFrame`` or flux-calibrated ``spCFrame`` files using the ``--frame`` or ``--cframe`` options, respectively.  These files contain a half plate of spectra for a single band (blue/red) and exposure.  Finally, you can load the ``spPlate`` files containing combined spectra for a whole plate using the ``--platefile`` option. See the :doc:`/overview` for details.

The ``--verbose`` option displays a progress bar showing the fraction of files already locally available. Any files that were previously fetched will not be downloaded again so it is safe and efficient to run ``bossfetch`` for overlapping lists of observations.  Note that the progress bar may appear to update unevenly if some files are already mirrored and others need to be downloaded.

Each data file download is streamed to a temporary files with ``.downloading`` appended to their name then renamed to remove this extension after the download completes normally. If a download is interrupted or fails for some reason, the partially downloaded file will remain in the local mirror.  Re-running a ``bossfetch`` command will automatically re-download any partially downloaded file.

By default, downloading is split between two parallel subprocesses but you can change this with the
``--nproc`` option.  For downloading "lite" files, using more than 2 subprocesses will probably not
improve the overall performance.

If you want to transfer large amounts of files, you should consider using `globus <https://www.globus.org>`_. To prepare a `globus` bulk data transfer file list, use the `--globus` option to specify the remote/local endpoint pair `remote#endpoint:local#endpoint`. Note that the `--save` option must also be used to specify an output filename. SDSS endpoints are documented at `here <http://www.sdss.org/dr12/data_access/bulk/>`_.

For example, to transfer files from `lbnl#sdss3` to `local#endpoint`::

    bossfetch qso.dat --globus lbnl#sdss3:username#endpoint --save globus-xfer.dat
    ssh username@cli.globusonline.org transfer -s 1 < globus-xfer.dat

.. _bossplot:

bossplot
--------

Plot the spectrum of a single BOSS observation, identified by its PLATE, MJD of the observation, and the FIBER that was assigned to the target whose spectrum you want to plot. For example (these are the defaults if you omit any parameters)::

    bossplot --plate 6641 --mjd 56383 --fiber 30

This should open a new window containing the plot that you will need to close in order to exit the program.  To also save your plot, add the ``--save-plot`` option with a filename that has a standard graphics format extension (pdf,png,...).  If you omit the filename, ``--save-plot`` uses the name ``bossplot-{plate}-{mjd}-{fiber}.png``. To save plots directly without displaying them, also use the ``--no-display`` option.

You can also save the data shown in a plot using ``--save-data`` with an optional filename (the default is ``bossplot-{plate}-{mjd}-{fiber}.dat``).  Data is saved using the `ascii.basic <http://docs.astropy.org/en/latest/api/astropy.io.ascii.Basic.html#astropy.io.ascii.Basic>`_ format and only wavelengths with valid data are included in the output.

Use ``--wlen-range [MIN:MAX]`` to specify a wavelength range over which to plot (x-axis), overriding the default, auto-detected range.  Similarly, ``--flux-range [MIN:MAX]`` and ``--wdisp-range [MIN:MAX]`` work for the flux (left y-axis) and dispersion (right y-axis).  MIN and MAX can be either blank (which means use the default value), an absolute value (1000), or a percentage (10%), and percentages and absolute values may be mixed.  Working examples::

    --wlen-range [:7500]
    --wlen-range [10%:90%]
    --wlen-range [10%:8000]

Note that a percentage value between 0-100% is interpreted as a percentile for vertical (flux, wdisp) axes. In all other cases, percentage values specify a limit value equal to a fraction of the full range ``[lo:hi]``::

    limit = lo + fraction*(hi - lo)

and can be < 0% or >100% to include padding. Another visual option ``--scatter`` will give a scatter plot of the flux rather than the flux 1-sigma error band.

Plots include a label ``PLATE-MJD-FIBER`` by default (or ``PLATE-MJD-FIBER-EXPID`` for a single exposure).  Add the option ``--label-pos <VALIGN>-<HALIGN>`` option to change its position, with ``<VALIGN> = top, center, bottom`` and ``<HALIGN> = left, center, right``.  Use ``--label-pos none`` to remove the label.  Use ``--no-grid`` to remove the default wavelength grid lines.

Several options are available to see data beyond just object flux.  Use ``--show-sky`` to show the subtracted sky (modeled) flux, ``--add-sky`` to show the total of object flux and modeled sky flux, ``--show-mask`` to show grayed regions where data has been masked out because it is deemed invalid, and ``--show-dispersion`` to show wavelength dispersion.

You will sometimes want to see data that would normally be masked as invalid. To include pixels with a particular `mask bit <http://www.sdss3.org/dr10/algorithms/bitmask_sppixmask.php>`_ set, use the ``--allow-mask`` option, e.g.::

    bossplot --allow-mask 'BRIGHTSKY|SCATTEREDLIGHT'

Note that multiple flags can be combined using the logical-or symbol ``|``, but this requires quoting as shown above. To show all data, including any invalid pixels, use the ``--show-invalid`` option.

The ``bossplot`` command will automatically download the appropriate data file if necessary.  This is 'conservative':  if an existing local file can be used to satisfy a request, no new files will be downloaded.

Spectra can be plotted from different data files. By default the spec-lite data file is used for a coadd or the spec file for an individual exposure.  Use the ``--frame`` or ``--cframe`` options to plot a single-exposure spectrum from a plate ``spFrame`` file or its flux-calibrated equivalent ``spCFrame`` file.  Use the ``--platefile`` option to plot the combined spectrum from an ``spPlate`` file. See the :doc:`/overview` for details.

To plot a single exposure, use the ``--exposure`` option to specify the sequence number (0,1,...) of the desired exposure. You can also set the ``--band`` option either ``blue`` or ``red`` to plot a single camera's data, or ``both`` to superimpose the overlapping data from both cameras.  Note that when displaying data from a co-added data product (spec, speclite, spPlate), the exposure sequence number only indexes exposures that were actually used in the final co-added spectrum.  However, the spFrame and spCFrame data products include all exposures used as input to the co-add (based on a :class:`bossdata.plate.Plan`) so, in cases where not all exposures are used, the ``--exposure`` option indexes a larger list of science exposures. Use the ``--verbose`` option to display information about the available exposures in either case.

This script uses the `matplotlib <http://matplotlib.org>`_ python library, which is not required for the ``bossdata`` package and therefore not automatically installed, but is included in scientific python distributions like `anaconda <https://store.continuum.io/cshop/anaconda/>`__.

.. _bossraw:

bossraw
-------

Assemble the raw CDD exposures used from one camera for a coadd into a single multi-extension FITS file.  For example::

    bossraw --plate 6641 --mjd 56383 --camera b1 --verbose

The MJD argument is optional in the common case that there is only one possible value.  The output will be saved to a file ``{plate}-{mjd}-{camera}.fits`` by default, or you can specify a file name with the ``--save`` argument.

The saved raw data is bias subtracted by default (or use ``--no-bias-subtraction``) and consists of only the :meth:`data regions <bossdata.raw.RawImageFile.get_amplifier_region>` of each CCD quadrant. See :meth:`bossdata.raw.RawImageFile.get_data` for details on the remaining optional arguments.

See the `sdR datamodel <http://data.sdss3.org/datamodel/files/BOSS_SPECTRO_DATA/MJD/sdR.html>`__ for more information about raw data.

To view the images stored in the output file, open it in `DS9 <http://ds9.si.edu/site/Home.html>`__ using the `File / Open As / Multiple Extension Frames...` menu item.
