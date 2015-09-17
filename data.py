#!/usr/bin/python

""" Load the Gaia-ESO Survey iDR4 node abundances. """

import numpy as np
from glob import glob
from collections import Counter

from astropy.table import Table

def get_elements(extension):
    """
    Return the elemental abundance columns available in the extension.
    """

    return { column: comment for _, column, comment in extension.header.cards \
        if "Abundance" in comment }
    


def retrieve_table(database, query, values=None, prefix=True):

    cursor = database.cursor()
    cursor.execute(query, values)
    names = [_[0] for _ in cursor.description]
    results = cursor.fetchall()
    if len(results) > 0:
        # Do we have multiple columns? If so, prefix them.
        duplicate_names = [k for k, v in Counter(names).items() if v > 1]
        if len(duplicate_names) > 0 and prefix:
            prefixes = []
            counted = []
            for name in names:
                if name in duplicate_names:
                    prefixes.append(counted.count(name))
                    counted.append(name)
                else:
                    prefixes.append(-1)
            names = [[n, "{0}.{1}".format(p, n)][p>=0] for p, n in zip(prefixes, names)]

        t = Table(rows=results, names=names)
        cursor.close()
        return t
    return None


def retrieve_column(database, query, values=None, asarray=True):

    cursor = database.cursor()
    cursor.execute(query, values)
    names = [_[0] for _ in cursor.description]
    results = cursor.fetchall()
    cursor.close()
    return results if not asarray else np.array(results).flatten()