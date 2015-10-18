============
Installation
============

To install, use the command line::

    % pip install bossdata

To upgrade to the latest version::

    % pip install bossdata --upgrade

Requirements
------------

The following additional pacakges are used by bossdata and will be installed automatically by pip, if necessary:

 * requests
 * progressbar
 * astropy
 * fitsio
 * numpy
 * pydl

Note that some operations run much slower (but still correctly) with numpy versions 10.0.0 and 10.0.1
so these should be avoided if possible. See `here <https://github.com/numpy/numpy/issues/6467>`__
for details.  To determine which version of numpy you are using::

    import numpy
    print numpy.version.version

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
