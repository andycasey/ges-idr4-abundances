#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Calculate and apply biases. """

__author__ = 'Andy Casey <arc@ast.cam.ac.uk>'

import logging
import numpy as np

import data
from utils import solar_abundance

logger = logging.getLogger('ges')


def _biases(measurements, column, e_column, zero_point):

    # Group by node, wavelength.
    biases = { node: {} for node in set(measurements["node"]) }
    measurements = measurements.group_by(["node", "wavelength"])
    for group in measurements.groups:
        node, wavelength = group["node"][0], group["wavelength"][0]

        # TODO weighted average?
        # TODO with covariance??
        mean = np.nanmean(group["abundance"])
        sigma = np.nanstd(group["abundance"])
        N = np.isfinite(group["abundance"]).sum()

        biases[node][wavelength] = (mean - zero_point, sigma, N)

    return biases


def solar_abundance_biases(database, element, ion, **kwargs):
    '''
    Calculate the abundance bias for each wavelength for each node, with respect
    to the Solar abundance.

    :param database:
        A PostgreSQL database connection.

    :type database:
        :class:`~psycopg2.connection`

    :param element:
        The element name to homogenise.

    :type element:
        str

    :param ion:
        The ionisation stage of the element to homogenise (1 = neutral).

    :type ion:
        int
    '''

    measurements = data.retrieve_table(database,
        """SELECT * FROM line_abundances WHERE element = '{0}' AND ion = '{1}'
        AND cname LIKE 'sssss%'""".format(element, ion))

    decimals = kwargs.pop("round_wavelengths", 1)
    if decimals >= 0:
        measurements["wavelength"] \
            = np.round(measurements["wavelength"], decimals)

    column = kwargs.pop("column", "abundance")
    e_column = kwargs.pop("e_column", "e_abundance")

    return _biases(measurements, column, e_column, solar_abundance(element))



def apply_solar_biases(database, element, ion):
    """

    """
    raise NotImplementedError



if __name__ == "__main__":

    import os
    kwds = {
        "host": "/tmp/",
        "dbname": "arc",
        "user": "arc",
        "password": os.environ.get('PSQL_PW', None)
    }


    import psycopg2 as pg
    database = pg.connect(**kwds)

    foo = solar_abundance_biases(database, "Si", 2)


