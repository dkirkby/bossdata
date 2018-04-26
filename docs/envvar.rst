Configuration
=============

You will normally want to establish your local configuration and specify which remote data you want to work with using some environment variables:

* ``BOSS_LOCAL_ROOT``: The top-level directory where all downloaded data files will be locally mirrored. Make sure there is enough space here for the files you plan to use locally. You might want to exclude this directory from your backups since it can get large and is already backed up remotely.
* ``BOSS_DATA_URL``: The top-level URL for downloading data, possibly including account information for accessing proprietary data.
* ``BOSS_SAS_PATH``: The top-level path of the data you want to work with, which will normally begin with "/sas".
* ``BOSS_REDUX_VERSION``: The pipeline reconstruction version that you want to work with.

If any of these variables is not specified, defaults appropriate for access the public `Data Release 14 <http://dr14.sdss.org>`_ will be used and any downloaded data will be saved to a temporary local directory. At a minimum, you should normally specify a permanent location for storing local data by setting the ``BOSS_LOCAL_ROOT`` environment variable.

The default settings of the other environment variables are equivalent to (in bash)::

    export BOSS_DATA_URL=http://dr14.sdss.org
    export BOSS_SAS_PATH=/sas/dr14/eboss
    export BOSS_REDUX_VERSION=v5_10_0

However these variables are set, the following unix shell command should always print a valid URL that displays a directory listing in any browser::

    echo $BOSS_DATA_URL/$BOSS_SAS_PATH/boss/spectro/redux/$BOSS_REDUX_VERSION/

You can optionally define one more environment variable ``BOSS_SPECLOG`` to locate a local checkout
of the ``speclog`` svn product.  This is only required if you need to access the full plug maps
(including non-science fibers) and prefer to use an environment variable instead of passing a
path argument.  See the :meth:`read_plug_map() <bossdata.raw.RawImageFile.read_plug_map>`
documentation for details.

If you are running on a system where the full dataset is available directly via the file system (e.g., at NERSC), use a ``$BOSS_DATA_URL`` that starts with the ``file://`` URI to indicate that data does not need to be transferred via network to your ``$BOSS_LOCAL_ROOT``.  In this case,
the local root will still be used for the sqlite files created by the :mod:`meta module <bossdata.meta>`. For details on using bossdata at NERSC, see :doc:`this guide </nersc>`.

The sections below describe how to access sources of data other than the default public DR12 release.

SDSS-III BOSS Data
------------------

To access the DR12 BOSS data release use::

    export BOSS_DATA_URL=http://dr12.sdss3.org
    export BOSS_SAS_PATH=/sas/dr12/boss
    export BOSS_REDUX_VERSION=v5_7_0

SEQUELS Data
------------

Quoting from `here <http://www.sdss.org/dr12/data_access/bulk/>`_:

    For BOSS, the main galaxy clustering survey is entirely contained in v5_7_0.
    After the main survey was finished, ancillary programs continued —
    these were processed as v5_7_2, which is the same code but a different
    processing version number to keep the datasets distinct.  The SEQUELS
    ancillary program has plates in both v5_7_0 and v5_7_2.

To access `SEQUELS <http://www.sdss.org/dr12/algorithms/ancillary/boss/sequels/>`_ data processed as v5_7_2, use::

    export BOSS_SAS_PATH=/sas/dr12/boss
    export BOSS_REDUX_VERSION=v5_7_2
    export BOSS_DATA_URL=http://dr12.sdss3.org

SDSS-I/II Spectra
-----------------

Some spectra from plates 0266 - 3006 are included in the public DR12 release and available under pipeline reduction versions 26, 103 and 104.  To access version 26, for example, use::

    export BOSS_SAS_PATH=/sas/dr12/sdss
    export BOSS_REDUX_VERSION=26
    export BOSS_DATA_URL=http://dr12.sdss3.org

eBOSS Proprietary Data
----------------------

Proprietary data from the `eBOSS survey <http://www.sdss.org/surveys/eboss/>`_ is password protected but still accessible via ``bossdata``.  Contact the authors for for details if you are an SDSS-IV collaborator.
