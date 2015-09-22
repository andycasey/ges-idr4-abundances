#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Calculate and apply biases. """

__author__ = 'Andy Casey <arc@ast.cam.ac.uk>'

import logging
import numpy as np
import scipy.optimize as op

import data
import utils

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

    return _biases(measurements, column, e_column, utils.solar_abundance(element))


def differential_abundance_biases(database, element, ion, scaled=False,
    ignore_flags=False, **kwargs):
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

    X, nodes, diff_data = utils.match_node_abundances(database, element, ion,
        scaled=scaled, ignore_flags=ignore_flags)

    # Calculate the full differential abundances.
    X_diff, indices = utils.calculate_differential_abundances(X, full_output=True)

    # Determine the differences to each node.
    diff_data["wavelength"] = np.round(diff_data["wavelength"], 1)
    wavelengths = set(diff_data["wavelength"])

    #node_masks = { n: (diff_data["node"] == n) for n in nodes }
    #wl_masks = { w: (diff_data["wavelength"] == w) for w in wavelengths }
    differential_bias = { n: { w: [] for w in wavelengths } for n in nodes }
    

    for wavelength in wavelengths:

        X_wl = X_diff[(diff_data["wavelength"] == wavelength), :]

        def differential_sigma(biases):
            print(biases)
            # Apply the biases on a per-column basis.
            X_offsets = np.zeros(X_wl.shape[1])
            for i, idx in enumerate(indices):

                # These are Node_0 - Node_1
                # We want to apply (Node_0 - offset_0) - (Node_1 - offset_1) so
                # the total offset is offset_1 - offset_0
                X_offsets[i] = biases[idx[1]] - biases[idx[0]]

            return np.nanstd(X_wl - X_offsets)

        result = op.fmin(differential_sigma, np.zeros(len(nodes)), disp=False)

        initial = differential_sigma(np.zeros(len(nodes)))
        final = differential_sigma(result)

        logger.info("Initial and final sigma: {0:.2f} {1:.2f}".format(initial,
            final))

        for node, offset in zip(nodes, result):
            differential_bias[node][wavelength] = (-offset, np.nan, -1)

    """
    for i, node in enumerate(nodes):
        for j, wavelength in enumerate(wavelengths):
            mask = node_masks[node] * wl_masks[wavelength]

            # Always make it Node - (all other nodes)
            values = np.hstack([(-1, +1)[i==idxs[0]] * X_diff[mask, k] \
                for k, idxs in enumerate(indices) if i in idxs]) 
            differential_bias[node][wavelength] = (np.nanmedian(values)/2.,
                np.nanstd(values), np.isfinite(values).sum())
    """
    return differential_bias




def apply_offset(database, element, ion, node, wavelength, offset,
    wavelength_tolerance=0.1):
    '''
    Apply an offset to the abundance measurements.
    '''

    if not np.isfinite(offset): return 0
    
    # Get the relevant line IDs
    ids = data.retrieve_column(database,
        """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = '{1}'
        AND node ILIKE '{2}%' AND wavelength > '{3}' AND wavelength < '{4}'"""\
        .format(element, ion, node, wavelength - wavelength_tolerance,
            wavelength + wavelength_tolerance), asarray=True)
    if len(ids) == 0: return 0

    # Apply the offset.
    rows = data.update(database,
        """UPDATE line_abundances
        SET scaled_abundance = abundance + {0} WHERE id IN %s""".format(offset),
        (tuple(ids), ))
    database.commit()
    return rows


if __name__ == "__main__":

    import os
    kwds = {
        "dbname": "arc",
    }


    import psycopg2 as pg
    database = pg.connect(**kwds)

    bar = differential_abundance_biases(database, "Si", 2)

    foo = solar_abundance_biases(database, "Si", 2)


