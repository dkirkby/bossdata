Using bossdata at NERSC
=======================

Use the following commands to install bossdata at NERSC, based on the DESI conda environment::

    # Use the DESI conda environment
    source /project/projectdirs/desi/software/desi_environment.sh

    # Setup a scratch area where additional packages can be installed
    cd $SCRATCH
    mkdir -p desi/lib/python3.6/site-packages desi/bin desi/code
    export PYTHONPATH=$SCRATCH/desi/lib/python3.6/site-packages:$PYTHONPATH
    export PATH=$SCRATCH/desi/bin:$PATH

    # Install bossdata into this scratch area
    cd $SCRATCH/desi/code
    git clone https://github.com/dkirkby/bossdata
    cd bossdata
    python setup.py develop --prefix $SCRATCH/desi

    # Create a directory for sqlite databases created by the meta module.
    mkdir -p $SCRATCH/bossdata

The commands are run once, for the initial bossdata installation. The following commands must be run each time you login to select this environment::

    # Use the DESI conda environment
    source /project/projectdirs/desi/software/desi_environment.sh

    # Use additional packages from the scratch area
    export PYTHONPATH=$SCRATCH/desi/lib/python3.6/site-packages:$PYTHONPATH
    export PATH=$SCRATCH/desi/bin:$PATH

    # Configure bossdata
    export BOSS_LOCAL_ROOT=$SCRATCH/bossdata
    export BOSS_DATA_URL=file:///global/projecta/projectdirs/sdss/www
    export BOSS_SAS_PATH=/sas/dr14/eboss
    export BOSS_REDUX_VERSION=v5_10_0

The environment variable settings above are for the public DR14 BOSS data release, but can be adjusted for other releases as described :doc:`here </envvar>`.  The key point is that your ``$BOSS_DATA_URL`` should start with ``file:///`` to indicate that all data files are available locally and so do not need to be downloaded via the network. In this case, ``$BOSS_LOCAL_ROOT`` is still used for the sqlite databases used by the :mod:`meta module <bossdata.meta>`.

To test your setup, try some of the :doc:`command line tools </scripts>`.

To update your version of bossdata to the latest master branch, use::

    cd $SCRATCH/desi/code/bossdata
    git pull
    python setup.py develop --prefix $SCRATCH/desi
