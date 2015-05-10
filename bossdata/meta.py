# -*- coding: utf-8 -*-

"""Support for querying the metadata associated with BOSS observations.
"""

from __future__ import division,print_function

import os.path
import gzip
import sqlite3

import astropy.table
import numpy as np

import bossdata.path
import bossdata.remote

def sql_create_table(table_name,recarray_dtype):
    """Prepare an SQL statement to create a database corresponding a numpy recarray data type.

    Assumes (but does not check) that columns named PLATE,MJD and FIBER are included and
    creates a composite primary index on these columns.

    Args:
        table_name(str): Name to give the new table.
        recarray_dtype: Numpy recarray data type that defines the columns to create.

    Returns:
        str: SQL statement that can be executed to create the database.

    Raises:
        ValueError: Cannot map data type to SQL.
    """
    columns = [ ]
    for (name,dtype) in recarray_dtype.descr:
        if np.issubdtype(dtype,np.integer):
            sql_type = 'INTEGER'
        elif np.issubdtype(dtype,np.float):
            sql_type = 'REAL'
        elif np.issubdtype(dtype,np.str):
            sql_type = 'TEXT'
        else:
            raise ValueError('Cannot map data type {} of {} to SQL.'.format(dtype,name))
        columns.append('`{name}` {type}'.format(name = name,type = sql_type))
    # Add a composite primary key on (plate,mjd,fiber).
    columns.append('PRIMARY KEY (PLATE,MJD,FIBER)')
    # Put the pieces together into the final SQL.
    return 'CREATE TABLE `{name}` ({columns})'.format(name = table_name,columns = ','.join(columns))

def create_meta_lite(sp_all_path,db_path,verbose = True):
    """Create the "lite" meta database from a locally mirrored spAll file.

    Args:
        sp_all_path(str): Absolute local path of the "lite" spAll file, which is expected to be
            a gzipped ASCII data file.
        db_path(str): Local path where the corresponding sqlite3 database will be written.
    """
    # Read the database into memory.
    if verbose:
        print('Reading...')
    with gzip.open(sp_all_path,mode = 'r') as f:
        table = astropy.table.Table.read(f,format = 'ascii')

    # Create a new database file.
    sql = sql_create_table('meta',table.dtype)
    connection = sqlite3.connect(db_path)
    cursor = connection.cursor()
    cursor.execute(sql)

    # Insert rows into the database.
    sql = 'INSERT INTO meta VALUES ({})'.format(','.join('?'*len(table.colnames)))
    if verbose:
        print('Writing...')
    for row in table:
        cursor.execute(sql,row)
    connection.commit()
    connection.close()

sql_type_map = {
    'INTEGER': np.integer,
    'REAL': float,
    'TEXT': str
}

class Database(object):
    """Initialize a searchable database of BOSS observation metadata.

    Args:
        finder(bossdata.path.Finder): Object used to find the names of BOSS data files.
        mirror(bossdata.remote.Manager): Object used to interact with the local mirror of BOSS data.
        lite(bool): Use the "lite" metadata format, which is considerably faster but only
            provides a subset of the most commonly accessed fields.
    """
    def __init__(self,finder,mirror,lite=True):

        # Get the local name of the metadata source file and the corresponding SQL database name.
        remote_path = finder.get_sp_all_path(lite = lite)
        local_path = mirror.local_path(remote_path)
        if lite:
            assert local_path.endswith('.dat.gz'),'Expected .dat.gz extension for {}.'.format(local_path)
            db_path = local_path.replace('.dat.gz','-lite.db')
        else:
            assert local_path.endswith('.fits'),'Expected .fits extention for {}.'.format(local_path)
            db_path = local_path.replace('.fits','.db')

        # Create the database if necessary.
        if not os.path.isfile(db_path):
            local_path = mirror.get(remote_path)
            if lite:
                create_meta_lite(local_path,db_path)
            else:
                raise RuntimeError('Creation of full meta db not supported yet.')

        # Connect to the database.
        self.connection = sqlite3.connect(db_path)
        self.cursor = self.connection.cursor()

        # Return TEXT values as ASCII strings when possible, instead of unicode.
        self.connection.text_factory = sqlite3.OptimizedUnicode

        # Look up and save the column definitions.
        self.cursor.execute('PRAGMA table_info(`meta`)')
        self.column_names = [ ]
        self.column_dtypes = [ ]
        for column_def in self.cursor:
            index,name,dtype = column_def[:3]
            self.column_names.append(str(name))
            self.column_dtypes.append(sql_type_map[dtype])

    def prepare_columns(self,column_names):
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
            return self.column_names,self.column_dtypes

        names,dtypes = [ ],[ ]
        for name in column_names.split(','):
            try:
                index = self.column_names.index(name)
            except ValueError:
                raise ValueError('Invalid column name {}.'.format(name))
            names.append(name)
            dtypes.append(self.column_dtypes[index])
        return names,dtypes

    def select_each(self,what = '*',where = None):
        """Iterate over the results of an SQL select query.

        This method is normally used as an iterator, e.g.

            for row in select(...):
                # each row is a tuple of values
                ...

        Since this method does not load all the results of a large query into memory, it is suitable
        for queries that are expected to return a large number of rows. For smaller queries, the
        :meth:`select_all` method might be more convenient.

        Args:
            what(str): Comma separated list of column names to return or '*' to return all columns.
            where(str): SQL selection clause or None for no filtering. Reserved column names such
                as PRIMARY must be `escaped` in this clause.
        """
        # Prepare the SQL select statement we will use.
        names,dtypes = self.prepare_columns(what)
        what = ','.join(['`{}`'.format(name) for name in names])
        sql = 'SELECT {} from meta'.format(what)
        if where:
            sql += ' WHERE {}'.format(where)

        # Execute the query and iterate over the results.
        self.cursor.execute(sql)
        for row in self.cursor:
            yield row

    def select_all(self,what = '*',where = None,max_rows = 100000):
        """Fetch all results of an SQL select query.

        Since this method loads all the results into memory, it is not suitable for queries that
        are expected to return a large number of rows.  Instead, use :meth:`select_each` for large
        queries.

        Args:
            what(str): Comma separated list of column names to return or '*' to return all columns.
            where(str): SQL selection clause or None for no filtering. Reserved column names such
                as PRIMARY must be `escaped` in this clause.
            max_rows(int): Maximum number of rows that will be returned.

        Returns:
            :class:`astropy.table.Table`: Table of results with column names matching those in
                the database, and column types inferred automatically.
        """
        # Prepare the SQL select statement we will use.
        names,dtypes = self.prepare_columns(what)
        what = ','.join(['`{}`'.format(name) for name in names])
        sql = 'SELECT {} from meta'.format(what)
        if where:
            sql += ' WHERE {}'.format(where)
        if max_rows:
            sql += ' LIMIT {:d}'.format(max_rows)

        # Execute the query and fetch all of the result into memory.
        self.cursor.execute(sql)
        rows = self.cursor.fetchall()

        # Return a table of the results, letting astropy.table.Table infer the columns types
        # from the data itself.
        return astropy.table.Table(rows = rows,names = names)
