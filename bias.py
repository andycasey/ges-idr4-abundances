#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Calculate and apply biases. """

__author__ = 'Andy Casey <arc@ast.cam.ac.uk>'

import logging
import numpy as np
import scipy.optimize as op

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


class AbundanceBiases(object):

    def __init__(self, release):
        self.release = release

    def solar(self, element, ion, **kwargs):
        """
        Calculate the abundance bias for each wavelength for each node, with
        respect to the Solar abundance.

        :param element:
            The element name to homogenise.

        :type element:
            str

        :param ion:
            The ionisation stage of the element to homogenise (1 = neutral).

        :type ion:
            int
        """

        measurements = self.release.retrieve_table(
            """SELECT * FROM line_abundances
            WHERE element = '{0}' AND ion = '{1}' AND cname LIKE 'sssss%'"""\
            .format(element, ion))

        decimals = kwargs.pop("round_wavelengths", 1)
        if decimals >= 0:
            measurements["wavelength"] \
                = np.round(measurements["wavelength"], decimals)

        column = kwargs.pop("column", "abundance")
        e_column = kwargs.pop("e_column", "e_abundance")

        return _biases(measurements, column, e_column, 
            utils.solar_abundance(element))


    def differential(self, element, ion, scaled=False, ignore_flags=False,
        **kwargs):
        """
        Calculate the differential abundance bias for each wavelength for each 
        node.

        :param element:
            The element name to homogenise.

        :type element:
            str

        :param ion:
            The ionisation stage of the element to homogenise (1 = neutral).

        :type ion:
            int
        """

        X, nodes, diff_data = utils.match_node_abundances(self.release._database,
            element, ion, scaled=scaled, ignore_flags=ignore_flags)

        # Calculate the full differential abundances.
        X_diff, indices = utils.calculate_differential_abundances(X,
            full_output=True)

        # Determine the differences to each node.
        diff_data["wavelength"] = np.round(diff_data["wavelength"], 1)
        wavelengths = set(diff_data["wavelength"])

        bias = { n: { w: (0, np.nan, -1) for w in wavelengths } for n in nodes }
        for wavelength in wavelengths:

            X_wl = X_diff[(diff_data["wavelength"] == wavelength), :]

            finite = { node: 0 for node in nodes }
            for k, (i, j) in enumerate(indices):
                value = np.isfinite(X_wl[:, k]).sum()
                finite[nodes[i]] += value
                finite[nodes[j]] += value

            finite_nodes = [node for node in nodes if finite[node] > 0]

            def differential_sigma(biases):
                # Apply the biases on a per-column basis.
                X_offsets = np.zeros(X_wl.shape[1])
                for i, idx in enumerate(indices):

                    # These are Node_0 - Node_1
                    # We want to apply (Node_0 - offset_0) - (Node_1 - offset_1)
                    # so the total offset is offset_1 - offset_0
                    # The biases.size is related to finite_nodes, not nodes.
                    try:
                        jdx0 = finite_nodes.index(nodes[idx[0]])
                        jdx1 = finite_nodes.index(nodes[idx[1]])
                    except ValueError:
                        continue
                    else:
                        X_offsets[i] = biases[jdx1] - biases[jdx0]

                return np.nanstd(X_wl - X_offsets)

            result = op.fmin(differential_sigma, np.zeros(len(finite_nodes)), 
                disp=False)
            
            initial = differential_sigma(np.zeros(len(finite_nodes)))
            final = differential_sigma(result)

            logger.info("Initial and final sigma: {0:.2f} {1:.2f}".format(
                initial, final))

            for node, offset in zip(finite_nodes, result):
                bias[node][wavelength] = (-offset, np.nan, -1)

        return bias


    def apply_offset(self, element, ion, node, wavelength, offset,
        wavelength_tolerance=0.1):
        """
        Apply an abundance offset to measurements of the given species (from
        the given line, by the given node). The offset is positively applied.

        :param element:
            The element name to homogenise.

        :type element:
            str

        :param ion:
            The ionisation stage of the element to homogenise (1 = neutral).

        :type ion:
            int

        :param node:
            The name of the node that the offset applies to.

        :type node:
            str

        :param wavelength:
            The wavelength of the line that the offset applies to.

        :type wavelength:
            float

        :param offset:
            The magnitude of the abundance offset to apply. This will be applied
            positively.

        :type offset:
            float

        :param wavelength_tolerance: [optional]
            The acceptable tolerance in the wavelength of the line.

        :type wavelength_tolerance:
            float

        :returns:
            The number of rows affected.
        """

        if not np.isfinite(offset): return 0
        
        # Get the relevant line IDs
        ids = self.release.retrieve_column("""SELECT id FROM line_abundances
            WHERE element = '{0}' AND ion = '{1}' AND node ILIKE '{2}%' 
            AND wavelength > '{3}' AND wavelength < '{4}'""".format(
                element, ion, node, wavelength - wavelength_tolerance,
                wavelength + wavelength_tolerance), asarray=True)
        if len(ids) == 0: return 0

        # Apply the offset.
        rows = self.release.update("""UPDATE line_abundances
            SET scaled_abundance = abundance + {0}
            WHERE id IN %s""".format(offset), (tuple(ids), ))
        self.release._database.commit()
        return rows