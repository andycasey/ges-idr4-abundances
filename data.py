#!/usr/bin/python

""" Load the Gaia-ESO Survey iDR4 node abundances. """

import logging
import numpy as np
from collections import Counter
from time import time

import psycopg2 as pg
from astropy.table import Table

# Initiate logging.
logging.basicConfig(level=logging.DEBUG,
    format='%(asctime)-15s %(name)s %(levelname)s %(message)s')
logger = logging.getLogger("ges")


def update(database, query, values=None, full_output=False, **kwargs):

    raise NotImplementedError

def execute(database, query, values=None, full_output=False, **kwargs):
    """
    Execute some SQL from the database.

    :param database:
        A PostgreSQL database connection.

    :type database:
        :class:`~psycopg2.connection`

    :param query:
        The SQL query to execute.

    :type query:
        str

    :param values: [optional]
        Values to use when formatting the SQL string.

    :type values:
        tuple or dict
    """
    
    t_init = time()
    try:    
        with database.cursor() as cursor:
            cursor.execute(query, values)
            results = cursor.fetchall()

    except pg.ProgrammingError:
        cursor.close()
        raise
    
    else:
        logger.info("Took {0:.0f} ms for SQL query {1}".format(
            1e3 * (time() - t_init), " ".join((query % values).split())))

    if full_output:
        names = tuple([column[0] for column in cursor.description])
        return (names, results, cursor.rowcount)

    return results


def retrieve_table(database, query, values=None, prefixes=True):
    """
    Retrieve a named table from a database.

    :param database:
        A PostgreSQL database connection.

    :type database:
        :class:`~psycopg2.connection`

    :param query:
        The SQL query to execute.

    :type query:
        str

    :param values: [optional]
        Values to use when formatting the SQL string.

    :type values:
        tuple or dict

    :param prefixes: [optional]
        Prefix duplicate column names with the given tuple.

    :type prefixes:
        tuple of str
    """

    names, rows, rowcount = execute(database, query, values, full_output=True)

    # TODO:
    if len(rows) == 0: return None

    counted_names = Counter(names)
    duplicates = [k for k, v in counted_names.items() if v > 1]
    if duplicates and prefixes:

        use_prefixes = map(str, range(max(counted_names.values()))) \
            if isinstance(prefixes, bool) else prefixes

        # Put the prefixes and names in the right order and format for joining.
        prefixes = [([], [use_prefixes[names[:i].count(n)]])[n in duplicates] \
            for i, n in enumerate(names)]
        names = [[n] for n in names]
        names = [".".join(p + n) for p, n in zip(prefixes, names)]

    return Table(rows=rows, names=names)


def retrieve_column(database, query, values=None, asarray=False):
    """
    Retrieve a single (unnamed) column from a database.

    :param database:
        A PostgreSQL database connection.

    :type database:
        :class:`~psycopg2.connection`

    :param query:
        The SQL query to execute.

    :type query:
        str

    :param values: [optional]
        Values to use when formatting the SQL string.

    :type values:
        tuple or dict

    :param asarray: [optional]
        Return the data as a numpy array.

    :type asarray:
        bool
    """

    rows = execute(database, query, values)
    return rows if not asarray else np.array(rows).flatten()