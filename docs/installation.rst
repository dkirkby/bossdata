============
Installation
============

To install, use the command line::

    % pip install bossdata

To upgrade to the latest version::

    % pip install bossdata --upgrade

Requirements
------------

The following additional packages are used by bossdata and will be installed automatically by pip, if necessary:

 * requests
 * progressbar2
 * astropy
 * fitsio
 * numpy
 * pydl

Numpy Performance Issue
^^^^^^^^^^^^^^^^^^^^^^^

Note that some operations run much slower (but still correctly) with numpy versions 1.10.0 and 1.10.1
so these should be avoided if possible. See `here <https://github.com/numpy/numpy/issues/6467>`__
for details.  To determine which version of numpy you are using::

    import numpy
    print numpy.version.version

The best solution is to use version 10.0.2 or later.  If this is not possible, revert to numpy 1.9.3 and astropy 1.0.4.  For example, with conda::

    conda install numpy=1.9.3
    conda install astropy=1.0.4

Optional Dependencies
---------------------

The following packages are optional and enable additional functionality.  They will not be
automatically installed by pip, but will be used when available.

 * matplotlib (used by the bossdata.plot module and bossplot script)

Quick Demonstration
-------------------

If you have `matplotlib <http://matplotlib.org>`_ installed, you can quickly test that everything is working with::

    bossplot

This should download a small data file for a single spectrum and plot the data in a window. Close the plot window to exit.  For more information on ``bossplot`` and other available command-line scripts, see :doc:`/scripts`.
