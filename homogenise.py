#!/usr/bin/python

""" Homogenise chemical abundances. """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import logging
import numpy as np
import scipy.optimize as op

import utils

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

        scaled = kwargs.get("scaled", True)

        # Remove any existing homogenised line or average abundances.
        self.release.execute("""DELETE FROM homogenised_line_abundances
            WHERE TRIM(element) = %s AND ion = %s""", (element, ion))
        self.release.execute("""DELETE FROM homogenised_abundances
            WHERE TRIM(element) = %s AND ion = %s""", (element, ion))

        # Drop index if it exists.
        self.release.execute(
            "DROP INDEX IF EXISTS homogenised_line_abundances_species_index")
        self.release.commit()

        # Get the unique wavelengths.
        wavelengths = sorted(set(self.release.retrieve_table(
            """SELECT DISTINCT ON (wavelength) wavelength FROM line_abundances
            WHERE TRIM(element) = %s AND ion = %s AND flags = 0
            ORDER BY wavelength ASC""", (element, ion))["wavelength"]))

        # Get the unique CNAMEs. We deal with repeat spectra (of the same CNAME)
        # later on in the code.
        cnames = self.release.retrieve_table(
            """SELECT DISTINCT ON (cname) cname FROM line_abundances
            WHERE TRIM(element) = %s AND ion = %s ORDER BY cname ASC""",
            (element, ion))["cname"]

        # For each wavelength, approximate the covariance matrix then homogenise
        # this wavelength for all cnames.
        column = "scaled_abundance" if scaled else "abundance"
        logger.debug("Homogenising {0} {1} using column {2}".format(element,
            ion, column))

        # In order to build the covariance matrix for this species, at some 
        # point we will need to know the variance for each line. We can estimate
        # this from the variance in the distribution of differential abundances.

        # This effectively tells us how well this line could have been measured
        # (after accounting for systematics between different nodes).
        Y, Y_nodes, Y_table = self.release._match_species_abundances(
            element, ion, scaled=scaled, include_flagged_lines=False)

        # Homogenise the line abundances for each wavelength separately.
        for i, wavelength in enumerate(sorted(set(wavelengths))):

            # Match all of the abundances for this given line, so that we can
            # use the array to calculate the correlation coefficients between
            # different nodes for this particular line.
            X, X_nodes = self.release._match_line_abundances(element, ion, 
                wavelength, column, ignore_gaps=True, include_limits=False,
                include_flagged_lines=False, **kwargs)

            # Get the measurement variance for this line.
            Z = utils.calculate_differential_abundances(
                Y[(Y_table["wavelength"] == wavelength)])

            line_variance = np.nanvar(np.abs(Z))

            if not np.isfinite(line_variance):
                # If the line variance is not finite, it means we do not have
                # differential abundances for this line. This is typically
                # because there were not enough nodes measuring this wavelength.

                # But we do know where this line sits with respect to the mean
                # abundance for this element. So we can still estimate the line
                # variance.
                assert Y.shape[1] > 1

                 # For each measured wavelength, what is the mean abundance for
                # the corresponding star, and what is the variance in that
                # distribution?
                matchers = []
                wl_mask = (Y_table["wavelength"] == wavelength)

                for stub in Y_table["spectrum_filename_stub"][wl_mask]:
                    stub_mask = (Y_table["spectrum_filename_stub"] == stub)
                    
                    value = Y[stub_mask*wl_mask].flatten()
                    node_mask = np.isfinite(value)
                    value = value[node_mask]

                    matchers.extend(Y[stub_mask * ~wl_mask, node_mask] - value)

                line_variance = np.nanvar(np.abs(matchers))

                if line_variance == 0:
                    # This line was only measured by one node in some stars, and
                    # in those stars there are no other measurements of this
                    # element, so we have no basis for the variance in this line
                    # This is a fringe case, and we will just have to do 
                    # something reasonable:
                    line_variance = kwargs.get("default_variance", 0.1**2)
                    logger.warn("Using default variance of {0:.2f} for {1} {2}"\
                        " line at {3}".format(line_variance, element, ion, 
                            wavelength))

            assert np.isfinite(line_variance) and line_variance > 0

            # For each CNAME / FILENAME, homogenise this line.
            for j, cname in enumerate(cnames):
                
                # The line_abundances function will need the element, ion,
                # wavelength, the variance in the line measurement, and the
                # correlation coefficients between nodes (or the matrix to
                # produce them), and the cname to know where to put things.
                result = self.line_abundances(cname, element, ion, wavelength,
                    line_variance, X, X_nodes)

        # Need to commit before we can do the averaged results per star.
        self.release.commit()

        # Create an index to speed things up.
        # To prevent parallel problems, first check that the index has not been
        # created by a parallel homogenisation script.
        try:    
            self.release.execute("""CREATE INDEX
                homogenised_line_abundances_species_index
                ON homogenised_line_abundances (cname, element, ion)""")
            self.release.commit()

        except:
            self.release.execute("rollback")

        # To homogenise the spectrum abundances, we will need the correlation
        # coefficients between each line.

        # Match the homogenised line abundances on a per-star basis.
        Q, Q_wavelengths = self.release._match_homogenised_line_abundances(
            element, ion, ignore_gaps=False, include_limits=False)
        Q_rho = np.atleast_2d(np.ma.corrcoef(Q))

        for j, cname in enumerate(cnames):
            self.spectrum_abundances(element, ion, cname, rho=Q_rho,
                rho_wavelengths=Q_wavelengths, **kwargs)

        self.release.commit()

        # TODO what should we return?
        return None


    def line_abundances(self, cname, element, ion, wavelength, line_variance,
        X=None, X_nodes=None, **kwargs):
        """
        Homogenise the line abundances for a given element, ion, star and 
        wavelength.

        :param cname:
            The CNAME of the star to perform line abundance homogenisation for.

        :type cname:
            str

        :param element:
            The element name to homogenise.

        :type element:
            str

        :param ion:
            The ionisation stage of the element to homogenise (1 = neutral).

        :type ion:
            int

        :param wavelength:
            The wavelength of the line to homogenise abundances for.

        :type wavelength:
            float

        :param line_variance:
            The variance ($\sigma^2$) for this given line. This is estimated
            from the sample covariance for this line, from all the nodes.

        :type line_variance:
            float

        :param X: [optional]
            A matrix containing the common-pairs of line abundances for this
            given line. The shape should be (N_nodes, N_common_stars). This
            array will be used to calculate the correlation coefficients between
            different nodes, for this particular line.

        :type X:
            :class numpy.array:

        :param X_nodes: [optional]
            The corresponding node names for the X array. This is required if
            X is given.

        :type X_nodes:
            tuple of str
        """

        if X is not None:
            assert X_nodes is not None
            assert len(X_nodes) == X.shape[0]

        # Get all valid line abundances for this element, ion, cname.
        if self.release.config.wavelength_tolerance > 0:   
            measurements = self.release.retrieve_table("""SELECT * 
                FROM line_abundances
                WHERE TRIM(element) = %s AND ion = %s AND cname = %s 
                AND flags = 0 AND abundance <> 'NaN'
                AND wavelength >= %s AND wavelength <= %s ORDER BY node ASC""",
                (element, ion, cname,
                    wavelength - self.release.config.wavelength_tolerance,
                    wavelength + self.release.config.wavelength_tolerance))
        else:
            measurements = self.release.retrieve_table("""SELECT *
                FROM line_abundances WHERE TRIM(element) = %s AND ion = %s 
                AND cname = %s AND flags = 0 AND abundance <> 'NaN'
                AND wavelength = %s ORDER BY node ASC""",
                (element, ion, cname, wavelength))

        if measurements is None:
            logger.warn("No valid {0} {1} measurements at {2} for CNAME = {3}"\
                .format(element, ion, wavelength, cname))
            return None
        

        rho = np.atleast_2d(np.corrcoef(X))
        rho_nodes = X_nodes


        # Account for repeat spectra for a single star.
        results = {}
        measurements = measurements.group_by(["spectrum_filename_stub"])
        for i, group in enumerate(measurements.groups):

            mu, sigma, N_nodes, is_limit \
                = _homogenise_spectrum_line_abundances(group, line_variance,
                    rho, rho_nodes, **kwargs)

            stub = group["spectrum_filename_stub"][0]
            results[stub] = self.insert_line_abundance(element, ion, cname,
                stub, wavelength, mu, sigma, N_nodes, is_limit)


            
        return results


    def spectrum_abundances(self, element, ion, cname, rho=None,
        rho_wavelengths=None, **kwargs):
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
            """SELECT * FROM homogenised_line_abundances WHERE TRIM(element) = %s
            AND ion = %s AND cname = %s""", (element, ion, cname))
        if line_abundances is None: return
        line_abundances = line_abundances.group_by(["spectrum_filename_stub"])

        # TODO: Sigma clipping.


        for i, group in enumerate(line_abundances.groups):

            if np.all(group["upper_abundance"] == 1):
                # TODO: this assumes limits come from elements with only one
                # line
                self.insert_spectrum_abundance(element, ion,
                    group["cname"][0], group["spectrum_filename_stub"][0],
                    group["abundance"][0], group["e_abundance"][0],
                    len(group), 1, True)

            else:   
                N_lines = len(group)

                # Build the covariance matrix.
                cov = np.ones((N_lines, N_lines))
                if np.all(np.isfinite(rho)):
                    for i, wavelength_i in enumerate(group["wavelength"]):
                        for j, wavelength_j in enumerate(group["wavelength"]):
                            if i != j:
                                i_index = rho_wavelengths.index(wavelength_i)
                                j_index = rho_wavelengths.index(wavelength_j)
                                cov[i, j] *= rho[i_index, j_index]
                else:
                    logger.warn("No correlation coefficients for {0} {1}".format(
                        element, ion))

                cov[~np.isfinite(cov)] = 0

                # Build the covariance matrix.
                # Variances come from the individual homogenised lines.
                sigmas = group["e_abundance"]
                cov *= (np.repeat(sigmas, N_lines) \
                    * np.tile(sigmas, N_lines)).reshape((N_lines, -1))

                def unbiased_variance(w, cov):
                    return np.sum(np.repeat(w, w.size) * np.tile(w, w.size) \
                        * cov.flatten())

                # Calculate the optimal weights that will produce the lowest
                # variance. TODO There is an analytic solution to this...
                def min_variance(weights):
                    # All weights must sum to one.
                    weights = weights.copy()
                    final_value = 1. - np.sum(weights)
                    weights = np.append(weights, final_value)
                    if np.any(weights < 0):
                        return np.inf

                    return unbiased_variance(weights, cov)

                w = op.fmin(min_variance, (np.ones(N_lines)/N_lines)[:-1],
                    disp=False)
                w = np.append(w, 1 - np.sum(w))

                # Calculate the homogenised values using the optimal weights.
                x = np.array(group["abundance"])
                mu = np.sum(w * x)
                variance = unbiased_variance(w, cov)

                # Ensure things are sensible, otherwise the matrix of correlation
                # coefficients may not be correct.
                if variance < 0:
                    cov = np.ones((N_lines, N_lines)) \
                        * (np.repeat(sigmas, N_lines) \
                        * np.tile(sigmas, N_lines)).reshape((N_lines, -1))

                    w = 1.0/np.diag(cov)
                    w /= w.sum()
                    mu = np.sum(w * x)
                    variance = unbiased_variance(w, cov)
                    logger.warn("Negative homogenised variance for {0} / {1}. "\
                        "Ignoring the correlation coefficients provided.".format(
                            group["cname"][0],
                            group["spectrum_filename_stub"][0]))

                assert variance > 0
                assert np.all(np.isfinite([mu, variance]))

                result = self.insert_spectrum_abundance(element, ion,
                    group["cname"][0], group["spectrum_filename_stub"][0],
                    mu, np.sqrt(variance), max(group["num_measurements"]),
                    N_lines, False)

        # TODO: should we average the abundances from multiple spectra?
        return None


    def insert_line_abundance(self, element, ion, cname, spectrum_filename_stub,
        wavelength, mu, sigma, N_nodes, is_limit, node_bitmask=0):
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

        :param mu:
            The homogenised line abundance.

        :type mu:
            float

        :param sigma:
            The 67th percentile (1-sigma uncertainty) in the homogenised line
            abundance.

        :type sigma:
            float

        :param N_nodes:
            The number of node measurements used to produce this line abundance.

        :type N_nodes:
            int

        :param is_limit:
            Whether the supplied values are an upper limit or not.

        :type is_limit:
            bool
        """

        # Does this record already exist? If so remove it.
        # TODO
        logger.warn(
            "Not checking for existing entries in homogenised_line_abundances")

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
                "abundance": mu,
                "e_abundance": sigma,
                "upper_abundance": int(is_limit),
                "num_measurements": N_nodes,
                "node_bitmask": node_bitmask
            })[0][0]


    def insert_spectrum_abundance(self, element, ion, cname, 
        spectrum_filename_stub, mu, sigma, N_nodes, N_lines, is_limit, **kwargs):
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

        :param mu:
            The homogenised line abundance.

        :type mu:
            float

        :param sigma:
            The 67th percentile (1-sigma uncertainty) in the homogenised
            abundance.

        :type sigma:
            float

        :param N_nodes:
            The number of node measurements used to produce this abundance.

        :type N_nodes:
            int

        :param N_lines:
            The number of lines used to produce this abundance.

        :type N_lines:
            int

        :param is_limit:
            Whether the supplied values are an upper limit or not.

        :type is_limit:
            bool
        """

        # Does this record already exist? If so remove it.
        # TODO
        logger.warn(
            "Not checking for existing entries in homogenised_abundances")

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
                "abundance": mu,
                "e_abundance": sigma,
                "num_measurements": N_nodes,
                "num_lines": N_lines,
                "upper_abundance": int(is_limit),
            })[0][0]


def _homogenise_spectrum_line_abundances(measurements, line_variance, rho=None,
    rho_nodes=None, **kwargs):

    # Note that by default we use scaled_abundance, not abundance, such that we
    # can apply any offsets or scaling as necessary *before* the homogenisation
    # begins.

    column = "scaled_abundance" if kwargs.get("scaled", True) else "abundance"
    N_nodes = np.isfinite(measurements[column]).sum()

    # We're deep here. Mark your assumptions!
    assert N_nodes > 0
    # One measurement per node.
    if len(set(measurements["node"])) != N_nodes:
        logger.warn("The number of unique nodes ({0}) does not match the "\
            "expected number ({1}): {2}".format(len(set(measurements["node"])),
                N_nodes, ", ".join(set(measurements["node"]))))

        # Ignore the other ones.
        ok = np.ones(len(measurements), dtype=bool)
        done = []
        for i, node in enumerate(measurements["node"]):
            if node in done: ok[i] = False
            else:
                done.append(node)
        logger.debug("Ignoring {0} rows".format(len(measurements) - ok.sum()))

        measurements = measurements[ok]
        N_nodes = len(set(measurements["node"]))

    assert np.isfinite(line_variance)

    # If we have measurements and upper limits, we will be conservative and
    # take the highest upper limit available.
    if np.any(measurements["upper_abundance"] > 0):
        is_limit = measurements["upper_abundance"] == 1
        limits = measurements[column][is_limit]
        N_nodes = len(set(measurements["node"]))
        # TODO: This assumes limits only exist when one line is present.
        return (np.nanmax(limits), np.nanstd(limits), N_nodes, True)

    # TODO: Sigma clipping.

    # Build the covariance matrix.
    cov = np.ones((N_nodes, N_nodes))
    if np.all(np.isfinite(rho)):
        for i, node_i in enumerate(measurements["node"]):
            for j, node_j in enumerate(measurements["node"]):
                if i != j:
                    i_index = rho_nodes.index(node_i)
                    j_index = rho_nodes.index(node_j)
                    cov[i, j] *= rho[i_index, j_index]

    else:
        logger.warn("Non-finite correlation coefficients for {0} {1} line at "\
            "{2:.0f} -- ignoring covariance between {3} nodes".format(
                measurements["element"][0], measurements["ion"][0],
                measurements["wavelength"][0], N_nodes))

    # Check that no node uncertainties are larger than the line variance
    # given to us. If they are, use that.
    variances = np.atleast_2d(np.array([
        [line_variance] * N_nodes,
        measurements["e_abundance"].astype(float)**2]))
    sigmas = np.sqrt(np.nanmax(variances, axis=0))

    # Build the covariance matrix.
    cov *= (np.repeat(sigmas, N_nodes) \
        * np.tile(sigmas, N_nodes)).reshape((N_nodes, -1))

    def unbiased_variance(w, cov):
        return np.sum(np.repeat(w, w.size) * np.tile(w, w.size) * cov.flatten())

    # Calculate the optimal weights that will produce the lowest variance.
    # TODO There is an analytic solution to this...
    def min_variance(weights):
        # All weights must sum to one.
        weights = weights.copy()
        final_value = 1. - np.sum(weights)
        weights = np.append(weights, final_value)

        if np.any(weights < 0):
            return np.inf

        return unbiased_variance(weights, cov)

    w = op.fmin(min_variance, (np.ones(N_nodes)/N_nodes)[:-1], disp=False)
    w = np.append(w, 1 - np.sum(w))

    # Calculate the homogenised values using the optimal weights.
    x = np.array(measurements[column])
    mu = np.sum(w * x)
    sigma = np.sqrt(unbiased_variance(w, cov))

    assert sigma > 0
    assert np.all(np.isfinite([mu, sigma]))

    return (mu, sigma, N_nodes, False)

