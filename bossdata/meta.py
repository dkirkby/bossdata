# -*- coding: utf-8 -*-
# Licensed under a MIT style license - see LICENSE.rst

"""Support for querying the metadata associated with BOSS observations.
"""

from __future__ import division, print_function

import os
import os.path
import gzip
import sqlite3
import stat
import re

import numpy as np
import astropy.table
import fitsio

from progressbar import ProgressBar, Percentage, Bar

import bossdata.path
import bossdata.remote


def sql_create_table(table_name, recarray_dtype, renaming_rules={}, primary_key=None):
    """Prepare an SQL statement to create a database for a numpy structured array.

    Any columns in the structured array data type that are themselves arrays will be
    unrolled to a list of scalar columns with names `COLNAME_I` for element [i] of a 1D array
    and `COLNAME_I_J` for element [i,j] of a 2D array, etc, with indices I,J,... starting
    from zero.

    Args:
        table_name(str): Name to give the new table.
        recarray_dtype: Numpy structured array data type that defines the columns to create.
        renaming_rules(dict): Dictionary of rules for renaming columns. There are no explicit
            checks that these rules do not create duplicate column names or that all rules
            are applied.
        primary_key(str): Column name(s) to use as the primary key, after apply renaming rules.
            No index is created if this argument is None.

    Returns:
        tuple: Tuple (sql,num_cols) where sql is an executable SQL statement to create the
            database and num_cols is the number of columns created.

    Raises:
        ValueError: Cannot map data type to SQL.
    """
    columns = []
    for column_dtype in recarray_dtype.descr:
        name, dtype = column_dtype[:2]
        # Rename this column if requested.
        if name in renaming_rules:
            name = renaming_rules[name]
        # Map this column's data type to one of the sqlite3 types.
        if np.issubdtype(dtype, np.integer):
            sql_type = 'INTEGER'
        elif np.issubdtype(dtype, np.float):
            sql_type = 'REAL'
        elif np.issubdtype(dtype, np.str):
            sql_type = 'TEXT'
        else:
            raise ValueError('Cannot map data type {} of {} to SQL.'.format(dtype, name))
        if len(column_dtype) == 2:
            columns.append('`{name}` {type}'.format(name=name, type=sql_type))
        else:
            # Handle sub-array columns.
            array_shape = column_dtype[2]
            array_ndim = len(array_shape)
            array_size = np.prod(array_shape)
            indices = np.unravel_index(np.arange(array_size), array_shape)
            for i in range(array_size):
                element_name = name
                for j in range(array_ndim):
                    element_name += '_{:d}'.format(indices[j][i])
                columns.append('`{name}` {type}'.format(name=element_name, type=sql_type))
    num_cols = len(columns)

    # Add a composite primary key on (plate,mjd,fiber).
    if primary_key:
        columns.append('PRIMARY KEY {}'.format(primary_key))
    # Put the pieces together into the final SQL.
    sql = 'CREATE TABLE `{name}` ({columns})'.format(
        name=table_name, columns=','.join(columns))
    return sql, num_cols


def create_meta_lite(sp_all_path, db_path, verbose=True):
    """Create the "lite" meta database from a locally mirrored spAll file.

    The created database has a composite primary index on the (PLATE,MJD,FIBER) columns and
    the input columns MODELFLUX0..4 are renamed MODELFLUX_0..4 to be consistent with their
    names in the full database after sub-array un-rolling.

    The DR12 spAll lite file is ~115Mb and converts to a ~470Mb SQL database file.
    The conversion takes about 3 minutes on a laptop with sufficient memory (~4 Gb).
    During the conversion, the file being written has the extension `.building` appended,
    then this extension is removed (and the file is made read only) once the conversion
    successfully completes.  This means that if the conversion is interrupted for any
    reason, it will be restarted the next time this function is called and you are
    unlikely to end up with an invalid database file.

    Args:
        sp_all_path(str): Absolute local path of the "lite" spAll file, which is expected
            to be a gzipped ASCII data file.
        db_path(str): Local path where the corresponding sqlite3 database will be written.
    """

    '''
    Things that have been tried RE: speeding this up to no (or very little) effect:
    1.)  Removing the primary key from the table creation, doing all insets, and
        then running SQL:
            CREATE UNIQUE INDEX psuedo_primary_key_meta ON meta (PLATE,MJD,FIBER)
        at very end.
    2.)  PRAGMA synchronous = OFF
        PRAGMA journal_mode = OFF
        and variations
    3.)  Explicitly setting transactions (which fails with an error, as this is
        being done implicitly)
    '''

    # Read the database into memory.
    if verbose:
        print('Initializing the lite database...')

    # Pass-1 uses only the first 100 lines to automatically determine the column names and
    # types.
    lines = ''
    with gzip.open(sp_all_path) as f:
        for i in range(100):
            lines += f.readline()
    table = astropy.table.Table.read(lines, format='ascii.fixed_width_two_line', guess=False)

    # Create a new database file.
    rules = {}
    for i in range(5):
        rules['MODELFLUX{}'.format(i)] = 'MODELFLUX_{}'.format(i)
    sql, num_cols = sql_create_table(
        'meta', table.dtype, renaming_rules=rules, primary_key='(PLATE,MJD,FIBER)')
    connection = sqlite3.connect(db_path + '.building')
    cursor = connection.cursor()
    cursor.execute(sql)
    connection.commit()

    data = np.loadtxt(sp_all_path, skiprows=2, dtype=table.dtype)

    # Insert rows into the database.
    sql = 'INSERT INTO meta VALUES ({values})'.format(values=','.join('?' * num_cols))
    if verbose:
        progress_bar = ProgressBar(
            widgets=['Writing', ' ', Percentage(), Bar()], maxval=len(data)).start()
    for i, row in enumerate(data):
        cursor.execute(sql, row)
        if verbose:
            progress_bar.update(i + 1)
    connection.commit()
    connection.close()

    # Make the temporary file read only by anyone.
    os.chmod(db_path + '.building', stat.S_IREAD | stat.S_IRGRP | stat.S_IROTH)
    # Move the temporary file to its permanent location.
    os.rename(db_path + '.building', db_path)

    if verbose:
        progress_bar.finish()


def create_meta_full(catalog_path, db_path, verbose=True, primary_key='(PLATE,MJD,FIBER)'):
    """Create the "full" meta database from a locally mirrored catalog file.

    The created database renames FIBERID to FIBER and has a composite primary index on the
    (PLATE,MJD,FIBER) columns. Sub-array columns are also unrolled: see
    :func:`sql_create_table` for details. The conversion takes about 24 minutes on a laptop
    with sufficient memory (~4 Gb). During the conversion, the file being written has the
    extension `.building` appended, then this extension is removed (and the file is made read
    only) once the conversion successfully completes.  This means that if the conversion is
    interrupted for any reason, it will be restarted the next time this function is called
    and you are unlikely to end up with an invalid database file.

    Args:
        catalog_path(str): Absolute local path of the "full" catalog file, which is expected
            to be a FITS file.
        db_path(str): Local path where the corresponding sqlite3 database will be written.
    """

    if verbose:
        print('Initializing the full database...')

    chunk_size = 100000
    position = 0

    # Open the FITs file.
    with fitsio.FITS(catalog_path) as hdulist:
        # This just reads the headers but still takes 15-20 seconds.
        # The equivalent operation with astropy.io.fits takes 15-20 minutes!

        # table = hdulist[1].read(rows=[position,position+chunk_size])
        table = hdulist[1][position:(position + chunk_size)]

        # Create a new database file.
        sql, num_cols = sql_create_table(
            'meta', table.dtype, renaming_rules={'FIBERID': 'FIBER'},
            primary_key=primary_key)
        connection = sqlite3.connect(db_path + '.building')
        cursor = connection.cursor()
        cursor.execute(sql)

        sql = 'INSERT INTO meta VALUES ({})'.format(','.join('?' * num_cols))
        if verbose:
            maxcount = hdulist[1].get_nrows()
            progress_bar = ProgressBar(
                widgets=['Writing', ' ', Percentage(), Bar()], maxval=maxcount).start()
            update_var = 1

        while len(table) > 0:
            # Insert rows into the database.
            for row in table:
                # Unroll columns with sub-arrays into a flat list to match the flat SQL schema,
                # and convert numpy types to the native python types required by sqlite3.
                values = []
                for j, column_data in enumerate(row):
                    if column_data.dtype.kind == 'S':
                        values.append(column_data.rstrip())
                    elif isinstance(column_data, np.ndarray):
                        values.extend(column_data.flatten().tolist())
                    else:
                        values.append(column_data.item())
                cursor.execute(sql, values)
                if verbose:
                    progress_bar.update(update_var)
                    update_var += 1

            position += chunk_size
            end_position = position + chunk_size
            if end_position > maxcount:
                end_position = maxcount
                if position > maxcount:
                    position = maxcount
            table = hdulist[1][position:end_position]

        connection.commit()
        connection.close()

        # Make the temporary file read only by anyone.
        os.chmod(db_path + '.building', stat.S_IREAD | stat.S_IRGRP | stat.S_IROTH)
        # Move the temporary file to its permanent location.
        os.rename(db_path + '.building', db_path)

        if verbose:
            progress_bar.finish()

sql_type_map = {
    'INTEGER': np.integer,
    'REAL': float,
    'TEXT': str
}


class Database(object):
    """Initialize a searchable database of BOSS observation metadata.

    Args:
        finder(bossdata.path.Finder): Object used to find the names of BOSS data files. If not
            specified, the default Finder constructor is used.
        mirror(bossdata.remote.Manager): Object used to interact with the local mirror of BOSS
            data. If not specified, the default Manager constructor is used.
        lite(bool): Use the "lite" metadata format, which is considerably faster but only
            provides a subset of the most commonly accessed fields. Ignored if either
            quasar_catalog or platelist is True.
        quasar_catalog(bool): Initialize database using the BOSS quasar catalog instead of
            spAll.
        quasar_catalog_name(str): The name of the BOSS quasar catalog to use, or use the
            :attr:`default <bossdata.path.Finder.default_quasar_catalog_name>` if this is None.
        platelist(bool): Initialize the database use the platelist catalog instead of spAll.
    """
    def __init__(self, finder=None, mirror=None, lite=True, quasar_catalog=False,
                 quasar_catalog_name=None, platelist=False, verbose=False):

        if finder is None:
            finder = bossdata.path.Finder(verbose=verbose)
        if mirror is None:
            mirror = bossdata.remote.Manager(verbose=verbose)

        # Get the local name of the metadata source file and the corresponding SQL
        # database name.
        if quasar_catalog:
            remote_path = finder.get_quasar_catalog_path(quasar_catalog_name)
            local_path = mirror.local_path(remote_path)
            assert local_path.endswith('.fits'), 'Expected .fits extention for {}.'.format(
                local_path)
            db_path = local_path.replace('.fits', '.db')
            lite_db_used = False
            # Create the database if necessary.
            if not os.path.isfile(db_path):
                local_path = mirror.get(remote_path)
                create_meta_full(local_path, db_path)
        elif platelist:
            remote_path = finder.get_platelist_path()
            local_path = mirror.local_path(remote_path)
            assert local_path.endswith('.fits'), 'Expected .fits extention for {}.'.format(
                local_path)
            db_path = local_path.replace('.fits', '.db')
            lite_db_used = False
            # Create the database if necessary.
            if not os.path.isfile(db_path):
                local_path = mirror.get(remote_path)
                create_meta_full(local_path, db_path, primary_key='(PLATE,MJD)')
        else:
            # Pre-build all our paths, test for (and store) the existence of the DB files
            remote_paths = [finder.get_sp_all_path(lite=True),
                            finder.get_sp_all_path(lite=False)]
            local_paths = [mirror.local_path(path) for path in remote_paths]
            db_paths = [Database._db_path_helper(local_paths[0], lite=True),
                        Database._db_path_helper(local_paths[1], lite=False)]
            db_paths_exist = [os.path.isfile(path) for path in db_paths]

            db_path = None
            local_path = None
            lite_db_used = True

            # Create the database if necessary.
            # Could make this more compact now that all the logic is foiled out.
            # lite is [0], full is [1]
            if lite and db_paths_exist[0]:        # lite branch, and lite DB exists
                db_path = db_paths[0]
            elif not lite and db_paths_exist[1]:  # full branch, and full DB exists
                db_path = db_paths[1]
                lite_db_used = False
            elif lite and not db_paths_exist[0]:  # lite branch and lite DB NOT exists
                if db_paths_exist[1]:             # ...but full does
                    db_path = db_paths[1]
                    lite_db_used = False
                else:                             # Neither DB's exist, so get files, create DB
                    local_path = mirror.get(remote_paths)
                    if local_path == local_paths[0]:    # lite
                        db_path = db_paths[0]
                        create_meta_lite(local_path, db_path)
                    else:                               # full
                        db_path = db_paths[1]
                        lite_db_used = False
                        create_meta_full(local_path, db_path)
            else:                                   # full branch and full DB NOT exists
                db_path = db_paths[1]
                lite_db_used = False
                local_path = mirror.get(remote_paths[1])
                create_meta_full(local_path, db_path)

        # Connect to the database.
        self.connection = sqlite3.connect(db_path)
        self.cursor = self.connection.cursor()

        # Return TEXT values as ASCII strings when possible, instead of unicode.
        self.connection.text_factory = sqlite3.OptimizedUnicode

        # Look up and save the column definitions.
        self.cursor.execute('PRAGMA table_info(`meta`)')
        self.column_names = []
        self.column_dtypes = []
        for column_def in self.cursor:
            index, name, dtype = column_def[:3]
            self.column_names.append(str(name))
            self.column_dtypes.append(sql_type_map[dtype])

        # Look up and save the number of rows in the database.
        self.cursor.execute('SELECT COUNT(*) FROM meta')
        self.num_rows = self.cursor.fetchone()[0]
        self.lite = lite_db_used

    def prepare_columns(self, column_names):
        """Validate column names and lookup their types.

        Args:
            column_names(str): Comma-separated list of column names or the special value '*' to
                indicate all available columns.

        Returns:
            tuple: Tuple (names,dtypes) of lists of column names and corresponding numpy
                data types.  Use :meth:`zip` to convert the return value into a recarray
                dtype.

        Raises:
            ValueError: Invalid column name.
        """
        if column_names == '*':
            return self.column_names, self.column_dtypes

        names, dtypes = [], []
        for name in re.split('\s*,\s*', column_names.strip()):
            try:
                index = self.column_names.index(name)
            except ValueError:
                raise ValueError('Invalid column name {}.'.format(name))
            names.append(name)
            dtypes.append(self.column_dtypes[index])
        return names, dtypes

    def select_each(self, what='*', where=None):
        """Iterate over the results of an SQL select query.

        This method is normally used as an iterator, e.g.

            for row in select(...):
                # each row is a tuple of values
                ...

        Since this method does not load all the results of a large query into memory, it
        is suitable for queries that are expected to return a large number of rows. For
        smaller queries, the :meth:`select_all` method might be more convenient.

        Args:
            what(str): Comma separated list of column names to return or '*' to return
                all columns.
            where(str): SQL selection clause or None for no filtering. Reserved column
                names such as PRIMARY must be escaped with backticks in this clause.

        Raises:
            sqlite3.OperationalError: failed to execute query.
        """
        # Prepare the SQL select statement we will use.
        names, dtypes = self.prepare_columns(what)
        what = ','.join(['`{}`'.format(name) for name in names])
        sql = 'SELECT {} from meta'.format(what)
        if where:
            sql += ' WHERE {}'.format(where)

        # Execute the query and iterate over the results.
        self.cursor.execute(sql)
        for row in self.cursor:
            yield row

    def select_all(self, what='*', where=None, sort=None, max_rows=100000):
        """Fetch all results of an SQL select query.

        Since this method loads all the results into memory, it is not suitable for queries
        that are expected to return a large number of rows.  Instead, use :meth:`select_each`
        for large queries.

        Args:
            what(str): Comma separated list of column names to return or '*' to return all
                columns.
            where(str): SQL selection clause or None for no filtering. Reserved column names
                such as PRIMARY must be escaped with backticks in this clause.
            max_rows(int): Maximum number of rows that will be returned.

        Returns:
            :class:`astropy.table.Table`: Table of results with column names matching those in
                the database, and column types inferred automatically. Returns None if no rows
                are selected.

        Raises:
            RuntimeError: failed to execute query.
        """
        # Prepare the SQL select statement we will use.
        names, dtypes = self.prepare_columns(what)
        what = ','.join(['`{}`'.format(name) for name in names])
        sql = 'SELECT {} from meta'.format(what)
        if where:
            sql += ' WHERE {}'.format(where)
        if sort:
            sql += ' ORDER BY {}'.format(sort)
        if max_rows:
            sql += ' LIMIT {:d}'.format(max_rows)

        # Execute the query and fetch all of the result into memory.
        try:
            self.cursor.execute(sql)
            rows = self.cursor.fetchall()
        except sqlite3.OperationalError as e:
            raise RuntimeError('Database error: {}'.format(e.args[0]))

        if len(rows) == 0:
            return None

        # Return a table of the results, letting astropy.table.Table infer the columns types
        # from the data itself.
        return astropy.table.Table(rows=rows, names=names)

    @staticmethod
    def _db_path_helper(path=None, lite=True):
        if lite:
            assert path.endswith('.dat.gz'), 'Expected .dat.gz extension for {}.'.format(
                path)
            return path.replace('.dat.gz', '-lite.db')
        else:
            assert path.endswith('.fits'), 'Expected .fits extention for {}.'.format(
                path)
            return path.replace('.fits', '.db')


def get_plate_mjd_list(plate, finder=None, mirror=None):
    """Return the list of MJD values when a plate was observed.

    Uses a query of the :datamodel:`platelist <platelist>`, so this file will be
    automatically downloaded if necessary.  Only MJD values for which the observation
    data quality is marked "good" will be returned.

    Args:
        plate(int): Plate number.
        finder(bossdata.path.Finder): Object used to find the names of BOSS data files. If not
            specified, the default Finder constructor is used.
        mirror(bossdata.remote.Manager): Object used to interact with the local mirror of BOSS
            data. If not specified, the default Manager constructor is used.

    Returns:
        list: A list of MJD values when this plate was observed. The list will
            be empty if this plate has never been observed.
    """
    platelist = Database(finder=finder, mirror=mirror, platelist=True)
    results = platelist.select_all(where='PLATE={0:d} and PLATEQUALITY="good"'.format(plate),
                                   what='MJD')
    if results:
        return list(results['MJD'])
    else:
        return []
