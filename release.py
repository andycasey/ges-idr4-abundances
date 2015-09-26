#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" A release object for convenience """

__author__ = 'Andy Casey <arc@ast.cam.ac.uk>'

import logging
from collections import Counter, namedtuple
from time import time

import numpy as np
import psycopg2 as pg
from astropy.table import Table
from matplotlib.cm import Paired

from bias import AbundanceBiases
from flag import AbundanceFlags
from homogenise import AbundanceHomogenisation
from plotting import AbundancePlotting

logging.basicConfig(level=logging.DEBUG,
    format='%(asctime)-15s %(name)s %(levelname)s %(message)s')
logger = logging.getLogger('ges')


class DataRelease(object):

    def __init__(self, database, user=None, password=None, host=None, **kwargs):
        """
        Initiate the data release object.

        """
        
        self._database = pg.connect(database=database,
            user=user, password=password, host=host, **kwargs)

        # Things we'll need.
        self.biases = AbundanceBiases(self)
        self.flags = AbundanceFlags(self)
        self.homogenise = AbundanceHomogenisation(self)
        self.plot = AbundancePlotting(self)

        Config = namedtuple('Configuration', 
            'wavelength_tolerance, round_wavelengths')

        self.config = Config(wavelength_tolerance=0.1, round_wavelengths=1)


    @property
    def nodes(self):
        try:
            return self._nodes
        except AttributeError:
            self._nodes = tuple(self.retrieve_column("""SELECT DISTINCT ON(node)
                node FROM line_abundances ORDER BY node ASC""", asarray=True))
            return self._nodes


    @property
    def node_colors(self):
        try:
            return self._node_colors
        except AttributeError:
            N = len(self.nodes)
            indices = np.linspace(0, 1, N)
            self._node_colors \
                = { n: Paired(indices[i]) for i, n in enumerate(self.nodes) }
            return self._node_colors


    def _match_species_abundances(self, element, ion, additional_columns=None,
        scaled=False, include_flagged_lines=False):
        """
        Return an array of matched line abundances for the given species.
        """

        column = "scaled_abundance" if scaled else "abundance"

        flag_query = "" if include_flagged_lines else "AND flags = 0"
        if additional_columns is None:
            data = self.retrieve_table("""SELECT node, wavelength,
                spectrum_filename_stub, abundance, scaled_abundance 
                FROM line_abundances WHERE trim(element) = %s AND ion = %s
                {flag_query} ORDER BY node ASC""".format(flag_query=flag_query),
                (element, ion))

        else:
            data = self.retrieve_table("""SELECT * FROM line_abundances l
                JOIN (SELECT DISTINCT ON (cname) cname, {additional_columns} 
                    FROM node_results ORDER BY cname) n 
                ON (trim(l.element) = %s AND l.ion = %s
                    AND l.cname = n.cname {flag_query})""".format(
                additional_columns=", ".join(additional_columns),
                flag_query=flag_query), (element, ion))

        if data is None: return (None, None, None)

        nodes = sorted(set(data["node"]))
        wavelengths = sorted(set(data["wavelength"]))

        data = data.group_by(["wavelength", "spectrum_filename_stub"])
        assert len(nodes) >= np.diff(data.groups.indices).max()

        N_nodes, N_groups = len(nodes), len(data.groups)
        X = np.nan * np.ones((N_groups, N_nodes))
        for i, group in enumerate(data.groups):
            logger.debug("matching group {0}/{1}".format(i, N_groups))
            for entry in group:
                j = nodes.index(entry["node"])
                X[i, j] = entry[column]

        # Return any other data.
        data = data[data.groups.indices[:-1]]
        return (X, nodes, data)


    def _match_line_abundances(self, element, ion, wavelength, column,
        ignore_gaps=False, include_flagged_lines=False):

        tol = self.config.wavelength_tolerance
        measurements = self.retrieve_table(
            "SELECT node, spectrum_filename_stub, " + column + \
            """ FROM line_abundances WHERE trim(element) = %s
            AND ion = %s {flag_query} AND wavelength >= %s
            AND wavelength <= %s""".format(
                flag_query="AND flags = 0" if not include_flagged_lines else ""),
                (element, ion, wavelength - tol, wavelength + tol))
        measurements = measurements.group_by(["spectrum_filename_stub"])

        nodes = sorted(set(measurements["node"]))
        X = np.nan * np.ones((len(nodes), len(measurements.groups)))
        for j, group in enumerate(measurements.groups):
            for row in group:
                X[nodes.index(row["node"]), j] = row[column]

        row_valid = np.all if ignore_gaps else np.any
        X = X[:, row_valid(np.isfinite(X), axis=0)]
        return (np.ma.array(X, mask=~np.isfinite(X)), nodes)


    def commit(self):
        """ Commit to the database. """
        return self._database.commit()


    def update(self, query, values=None, full_output=False, **kwargs):
        """
        Update the database with a SQL query.

        :param query:
            The SQL query to execute.

        :type query:
            str

        :param values: [optional]
            Values to use when formatting the SQL string.

        :type values:
            tuple or dict
        """

        logger.debug("Running SQL update query: {}".format(query))
        names, results, cursor = self.execute(query, values, **kwargs)
        return (names, results, cursor) if full_output else cursor.rowcount
        

    def retrieve(self, query, values=None, full_output=False, **kwargs):
        """
        Retrieve some data from the database.

        :param query:
            The SQL query to execute.

        :type query:
            str

        :param values: [optional]
            Values to use when formatting the SQL string.

        :type values:
            tuple or dict
        """

        names, results, cursor = self.execute(query, values, fetch=True,
            **kwargs)
        return (names, results, cursor.rowcount) if full_output else results


    def execute(self, query, values=None, fetch=False, **kwargs):
        """
        Execute some SQL from the database.

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
            with self._database.cursor() as cursor:
                cursor.execute(query, values)
                if fetch: results = cursor.fetchall()
                else: results = None

        except pg.ProgrammingError:
            logger.exception("SQL query failed:")
            cursor.close()
            raise
        
        else:
            taken = 1e3 * (time() - t_init)
            try:
                logger.info("Took {0:.0f} ms for SQL query {1}".format(taken,
                    " ".join((query % values).split())))
            except (TypeError, ValueError):
                logger.info("Took {0:.0f} ms for SQL query {1} with values {2}"\
                    .format(taken, query, values))

        names = None if cursor.description is None \
            else tuple([column[0] for column in cursor.description])
        return (names, results, cursor)


    def retrieve_table(self, query, values=None, prefixes=True):
        """
        Retrieve a named table from a database.

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

        names, rows, rowcount = self.retrieve(query, values, full_output=True)

        # TODO:
        if len(rows) == 0:
            return None

        counted_names = Counter(names)
        duplicates = [k for k, v in counted_names.items() if v > 1]
        if duplicates and prefixes:

            use_prefixes = map(str, range(max(counted_names.values()))) \
                if isinstance(prefixes, bool) else prefixes

            # Put the prefixes and names in the right order & format for joining
            prefixes = [
                ([], [use_prefixes[names[:i].count(n)]])[n in duplicates] \
                for i, n in enumerate(names)]
            names = [[n] for n in names]
            names = [".".join(p + n) for p, n in zip(prefixes, names)]

        table = Table(rows=rows, names=names)
        # If wavelengths are in the table, and we should round them, round them.
        if "wavelength" in table.dtype.names \
        and self.config.round_wavelengths >= 0:
            table["wavelength"] = np.round(table["wavelength"],
                self.config.round_wavelengths)
        return table


    def retrieve_column(self, query, values=None, asarray=False):
        """
        Retrieve a single (unnamed) column from the database.

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

        rows = self.retrieve(query, values)
        return rows if not asarray else np.array(rows).flatten()
