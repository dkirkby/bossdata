Executable scripts
==================

bossquery
---------

Query the meta data for BOSS observations. For example::

    bossquery --what PLUG_RA,PLUG_DEC,Z --where 'OBJTYPE="QSO"' --save qso.dat

This command uses an sqlite3 database of metadata that will be created if necessary. The save option supports `many different output formats<http://astropy.readthedocs.org/en/latest/io/unified.html#built-in-table-readers-writers>`_ that are automatically selected based on the file extension.  In addition, this program automatically maps the `.dat` and `.txt` extensions to the `ascii` format.
