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

Quick Demonstration
-------------------

If you have `matplotlib <http://matplotlib.org>`_ installed, you can quickly test that everything is working with::

    bossplot

This should download a small data file for a single spectrum and plot the data in a window. Close the plot window to exit.  For more information on ``bossplot`` and other available command-line scripts, see :doc:`/scripts`.
