#!/usr/bin/python

""" Homogenise chemical abundances. """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import logging
import numpy as np

logger = logging.getLogger("ges")


class AbundanceHomogenisation(object):

    def __init__(self, release):
        self.release = release


    def species(self, element, ion, **kwargs):
        """
        Homogenise the line abundances for all stars for a given element and ion.

        :param element:
            The element name to homogenise.

        :type element:
            str

        :param ion:
            The ionisation stage of the element to homogenise (1 = neutral).

        :type ion:
            int
        """

        # Remove any existing homogenised line or average abundances.
        self.release.execute("""DELETE FROM homogenised_line_abundances
            WHERE trim(element) = %s AND ion = %s""", (element, ion))
        self.release.execute("""DELETE FROM homogenised_abundances
            WHERE trim(element) = %s AND ion = %s""", (element, ion))
        self.release.commit()

        # Get the unique wavelengths.
        wavelengths = self.release.retrieve_column(
            """SELECT DISTINCT ON (wavelength) wavelength FROM line_abundances
            WHERE trim(element) = %s AND ion = %s AND flags = 0
            ORDER BY wavelength ASC""", (element, ion), asarray=True)

        if self.release.config.round_wavelengths >= 0:
            wavelengths = \
                np.round(wavelengths, self.release.config.round_wavelengths)

        # Get the unique CNAMEs.
        cnames = self.release.retrieve_column(
            """SELECT DISTINCT ON (cname) cname FROM line_abundances
            WHERE trim(element) = %s AND ion = %s ORDER BY cname ASC""",
            (element, ion), asarray=True)

        # For each wavelength, approximate the covariance matrix then homogenise
        # this wavelength for all cnames.
        column = ("abundance", "scaled_abundance")[kwargs.get("scaled", True)]
        logger.debug("Homogenising {0} {1} using column {2}".format(element,
            ion, column))

        for i, wavelength in enumerate(set(wavelengths)):
            # Calculate the correlation coefficients for this line.
            X, nodes = self.release._match_line_abundances(element, ion, 
                wavelength, column, ignore_gaps=True, 
                include_flagged_lines=False, **kwargs)
            rho = np.atleast_2d(np.corrcoef(X))
            cov = np.atleast_2d(np.cov(X))

            # For each cname, homogenise this line.
            for j, cname in enumerate(cnames):
                self.line_abundances(element, ion, cname, wavelength, rho, cov,
                    nodes, **kwargs)

        # Need to commit before we can do the averaged results per star.
        self.release.commit()

        for j, cname in enumerate(cnames):
            self.spectrum_abundances(element, ion, cname, **kwargs)

        self.release.commit()

        # TODO what should we return?
        return None


    def line_abundances(self, element, ion, cname, wavelength, rho, cov, nodes,
        **kwargs):
        """
        Homogenise the line abundances for a given element, ion, star and 
        wavelength.

        :param element:
            The element name to homogenise.

        :type element:
            str

        :param ion:
            The ionisation stage of the element to homogenise (1 = neutral).

        :type ion:
            int

        :param cname:
            The CNAME of the star to perform line abundance homogenisation for.

        :type cname:
            str

        :param wavelength:
            The wavelength of the line to homogenise abundances for.

        :type wavelength:
            float

        :param rho:
            The node-to-node correlation coefficients for this species.

        :type rho:
            :class:`~np.array`
        """

        # Get all valid line abundances for this element, ion, cname.
        if self.release.config.wavelength_tolerance > 0:   
            measurements = self.release.retrieve_table("""SELECT * 
                FROM line_abundances
                WHERE trim(element) = %s AND ion = %s AND cname = %s 
                AND flags = 0 AND abundance <> 'NaN'
                AND wavelength >= %s AND wavelength <= %s ORDER BY node ASC""",
                (element, ion, cname,
                    wavelength - self.release.config.wavelength_tolerance,
                    wavelength + self.release.config.wavelength_tolerance))
        else:
            measurements = self.release.retrieve_table("""SELECT *
                FROM line_abundances WHERE trim(element) = %s AND ion = %s 
                AND cname = %s AND flags = 0 AND abundance <> 'NaN'
                AND wavelength = %s ORDER BY node ASC""",
                (element, ion, cname, wavelength))
        if measurements is None: return
        
        # Account for repeat spectra for a single star.
        results = {}
        measurements = measurements.group_by(["spectrum_filename_stub"])
        for i, group in enumerate(measurements.groups):

            mean, sigma, N = _homogenise_spectrum_line_abundances(
                group, rho, cov, nodes, **kwargs)

            stub = group["spectrum_filename_stub"][0]
            results[stub] = self.insert_line_abundance(element, ion, cname,
                stub, wavelength, mean, sigma, N)
            
        return results


    def spectrum_abundances(self, element, ion, cname, **kwargs):
        """
        Average the homogenised line abundances for a given star.

        :param element:
            The element name to homogenise.

        :type element:
            str

        :param ion:
            The ionisation stage of the element to homogenise (1 = neutral).

        :type ion:
            int

        :param cname:
            The CNAME of the star to perform line abundance homogenisation for.

        :type cname:
            str
        """
        
        # Get the data for this star/spectrum filename stub.
        line_abundances = self.release.retrieve_table(
            """SELECT * FROM homogenised_line_abundances WHERE trim(element) = %s
            AND ion = %s AND cname = %s""", (element, ion, cname))
        if line_abundances is None: return
        line_abundances = line_abundances.group_by(["spectrum_filename_stub"])

        clip_outliers = kwargs.pop("clip_outliers", 5)

        for i, group in enumerate(line_abundances.groups):

            """
            # Any really bad outliers?
            if clip_outliers and clip_outliers is not None:
                sigma = np.nanstd(group["abundance"])
                outlier = np.abs(group["abundance"] - np.nanmean(group["abundance"]))
                is_outlier = outlier > clip_outliers

                # Remove this from the group
                group = group[~is_outlier]
            """

            # Calculate a weighted mean and variance.
            # TODO: currently ignoring the covariance between lines.
            inv_variance = 1.0/(group["e_abundance"]**2)
            weights = inv_variance/np.sum(inv_variance)

            mu = np.sum(group["abundance"] * weights)
            sigma = np.sqrt(np.sum(weights**2 * group["e_abundance"]**2))

            # Update the line abundance.
            assert len(group) == len(set(group["wavelength"]))
            N_lines = len(set(group["wavelength"]))
            N_measurements = np.sum(np.isfinite(group["abundance"]))

            result = self.insert_spectrum_abundance(element, ion,
                group["cname"][0], group["spectrum_filename_stub"][0],
                mu, sigma, N_lines, N_measurements)
        
        # TODO: what should we return?
        # TODO: should we average the abundances from multiple spectra?
        return None


    def insert_line_abundance(self, element, ion, cname, 
        spectrum_filename_stub, wavelength, abundance, uncertainty,
        num_measurements, node_bitmask=0):
        """
        Update the homogenised line abundance for a given CNAME, species
        (element and ion), and wavelength.

        :param element:
            The element name to homogenise.

        :type element:
            str

        :param ion:
            The ionisation stage of the element to homogenise (1 = neutral).

        :type ion:
            int

        :param cname:
            The CNAME of the star to perform line abundance homogenisation for.

        :type cname:
            str

        :param wavelength:
            The wavelength of the line.

        :type wavelength:
            float

        :param abundance:
            The homogenised line abundance.

        :type abundance:
            float

        :param uncertainty:
            The uncertainty in the homogenised line abundance.

        :type uncertainty:
            float

        :param num_measurements:
            The number of measurements used to produce this line abundance.

        :type num_measurements:
            int
        """

        # Does this record already exist? If so remove it.
        # TODO
        print("warning not checking for existing entries")

        return self.release.retrieve("""INSERT INTO homogenised_line_abundances
            (cname, spectrum_filename_stub, wavelength, element, ion, abundance,
            e_abundance, upper_abundance, num_measurements, node_bitmask)
            VALUES (%(cname)s, %(spectrum_filename_stub)s, %(wavelength)s,
            %(element)s, %(ion)s, %(abundance)s, %(e_abundance)s,
            %(upper_abundance)s, %(num_measurements)s, %(node_bitmask)s)
            RETURNING id;""", {
                "spectrum_filename_stub": spectrum_filename_stub,
                "cname": cname,
                "wavelength": wavelength,
                "element": element,
                "ion": ion,
                "abundance": abundance,
                "e_abundance": uncertainty,
                "upper_abundance": 0, # TODO
                "num_measurements": num_measurements,
                "node_bitmask": node_bitmask
            })[0][0]


    def insert_spectrum_abundance(self, element, ion, cname, 
        spectrum_filename_stub, abundance, uncertainty, num_lines, 
        num_measurements, **kwargs):
        """
        Update the homogenised abundance for a species (element and ion) as 
        measured in a spectrum (CNAME and spectrum filename stub).

        :param element:
            The element name to homogenise.

        :type element:
            str

        :param ion:
            The ionisation stage of the element to homogenise (1 = neutral).

        :type ion:
            int

        :param cname:
            The CNAME of the star to perform line abundance homogenisation for.

        :type cname:
            str

        :param spectrum_filename_stub:
            The spectrum filename stub for the given abundance.

        :type spectrum_filename_stub:
            str

        :param abundance:
            The homogenised line abundance.

        :type abundance:
            float

        :param uncertainty:
            The uncertainty in the homogenised line abundance.

        :type uncertainty:
            float

        :param num_lines:
            The number of transitions used to produce this line abundance.

        :type num_lines:
            int

        :param num_measurements:
            The number of measurements (from all nodes and transitions) used to
            produce this line abundance.

        :type num_measurements:
            int
        """

        # Does this record already exist? If so remove it.
        # TODO
        print("warning not checking for existing entries")

        return self.release.retrieve("""INSERT INTO homogenised_abundances
            (cname, spectrum_filename_stub, element, ion, abundance,
            e_abundance, upper_abundance, num_lines, num_measurements)
            VALUES (%(cname)s, %(spectrum_filename_stub)s,
            %(element)s, %(ion)s, %(abundance)s, %(e_abundance)s,
            %(upper_abundance)s, %(num_lines)s, %(num_measurements)s)
            RETURNING id;""", {
                "spectrum_filename_stub": spectrum_filename_stub,
                "cname": cname,
                "element": element,
                "ion": ion,
                "abundance": abundance,
                "e_abundance": uncertainty,
                "num_measurements": num_measurements,
                "num_lines": num_lines,
                "upper_abundance": 0, # TODO
            })[0][0]


def _homogenise_spectrum_line_abundances(measurements, rho, cov, nodes,
    **kwargs):

    # Note that by default we use scaled_abundance, not abundance, such that we
    # can apply any offsets or scaling as necessary *before* the homogenisation
    # begins.

    e_column = "e_abundance"
    column = ("abundance", "scaled_abundance")[kwargs.get("scaled", True)]

    N_measurements = np.isfinite(measurements[column]).sum()

    if N_measurements == 0:
        return (np.nan, np.nan, -1)

    # Problems to address later:
    assert np.all(measurements["upper_abundance"] == 0)

    abundance = np.nan * np.ones(len(nodes))
    u_abundance = np.nan * np.ones(len(nodes))
    for i, node in enumerate(nodes):
        if node not in measurements["node"]: continue
        match = (measurements["node"] == node)
        abundance[i] = measurements[column][match][0]
        #u_abundance[i]  = measurements[e_column][match][0] # use NODE ABUNDANCE
        u_abundance[i] = (np.diag(cov)[i]**0.5)/np.sqrt(len(nodes))
        if np.isfinite(measurements[e_column][match][0]):
            u_abundance[i] = np.sqrt(u_abundance[i]**2 + measurements[e_column][match][0]**2)

    mask = np.isfinite(abundance)
    abundance = abundance[mask]
    u_abundance = u_abundance[mask]

    # Fill in missing abundance uncertainties with the mean of the others.
    #bad_values = ~np.isfinite(u_abundance) + (0 >= u_abundance)
    #u_abundance[bad_values] = np.nanmean(u_abundance[~bad_values])
    if not np.all(np.isfinite(u_abundance)):
        logger.warn("Setting 0.2 dex uncertainty for {0} {1} line at {2}".format(
            None, None, measurements["wavelength"][0]))
        u_abundance[:] = 0.2

    N = mask.sum()
    N_nodes = len(nodes)
    mask = np.tile(mask, N_nodes).reshape(N_nodes, N_nodes)
    mask = np.multiply(mask, mask.T)

    sigmas = np.tile(u_abundance, N).reshape(N, N)
    sigmas = np.multiply(sigmas, sigmas.T)
    cov = np.multiply(rho[mask].reshape(N, N), sigmas)

    # Calculate the weights
    inv_variance = 1.0/(u_abundance**2)
    weights = inv_variance/np.sum(inv_variance)

    # Determine the weighted average and uncertainty.
    abundance_mean = np.sum(weights * abundance)
    abundance_sigma = np.sqrt(np.sum(
        np.tile(weights, N) * np.repeat(weights, N) * cov.flatten()))

    if not np.isfinite(abundance_sigma):
        abundance_sigma = np.nanstd(abundance)

    assert np.isfinite(abundance_mean * abundance_sigma)

    # Place this information in the homogenised abundance table.
    return (abundance_mean, abundance_sigma, N)
