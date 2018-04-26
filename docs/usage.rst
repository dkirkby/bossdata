=========
API Usage
=========

To use the ``bossdata`` package in your own python projects, you will normally start with::

    import bossdata
    finder = bossdata.path.Finder()
    mirror = bossdata.remote.Manager()

This code will use the environment variables ``$BOSS_SAS_PATH``, ``$BOSS_REDUX_VERSION``, ``$BOSS_DATA_URL`` and ``$BOSS_LOCAL_ROOT`` to configure your access to SDSS data files (see :doc:`/envvar` for details.) The ``finder`` and ``mirror`` objects can be used together to access locally mirrored copies of BOSS data files. For example::

    remote_path = finder.get_spec_path(plate=4567, mjd=55589, fiber=88, lite=True)
    local_path = mirror.get(remote_path)

Refer to the :doc:`API documentation </api>` for details on using the :mod:`bossdata.path` and :mod:`bossdata.remote` modules.

Certain data files have a helper class for accessing their contents:
 * spec,spec-lite: :class:`bossdata.spec.SpecFile`
 * plate: :class:`bossdata.plate.PlateFile`
 * plan: :class:`bossdata.plate.Plan`
 * frame,cframe: :class:`bossdata.plate.FrameFile`

For example, to open the spec-lite file used in the example above, use::

    spec = bossdata.spec.SpecFile(local_path)

The pattern for accessing large metadata files is somewhat different, and handled by the :class:`bossdata.meta.Database` class.
