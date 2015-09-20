#!/usr/bin/python

""" Homogenise line abundances. """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import logging
import numpy as np
from time import time

import data

logger = logging.getLogger("ges")

def create_table(database):
    """
    Create a table for the homogenised line abundances.

    :param database:
        A PostgreSQL database connection.

    :type database:
        :class:`~psycopg2.connection`
    """

    # Superfluous column removed:
    #object char(21) not null,
    
    cursor = database.cursor()
    cursor.execute("""CREATE TABLE homogenised_line_abundances (
        cname char(16) not null,
        spectrum_filename_stub char (140) not null,
        wavelength real not null,
        element char(2) not null,
        ion integer not null,
        abundance real,
        e_abundance real,
        upper_abundance int default 0,
        num_measurements int default 0,
        node_bitmask int default 0
        );""")
    cursor.execute("""ALTER TABLE homogenised_line_abundances 
        ADD COLUMN id BIGSERIAL PRIMARY KEY;""")
    cursor.close()

    cursor = database.cursor()
    cursor.execute("""CREATE TABLE homogenised_abundances (
        cname char(16) not null,
        spectrum_filename_stub char(140) not null,
        element char(2) not null,
        ion integer not null,
        abundance real,
        e_abundance real,
        upper_abundance int default 0,
        num_lines int default 0,
        num_measurements int default 0);""")
    cursor.execute("""ALTER TABLE homogenised_abundances ADD COLUMN id
        BIGSERIAL PRIMARY KEY;""")
    cursor.close()

    database.commit()
    return True


def drop_table(database):
    """
    Drop the table used for homogenised line abundances.

    :param database:
        A PostgreSQL database connection.

    :type database:
        :class:`~psycopg2.connection`
    """

    cursor = database.cursor()
    cursor.execute("DROP TABLE IF EXISTS homogenised_line_abundances;")
    cursor.execute("DROP TABLE IF EXISTS homogenised_abundances;")
    cursor.close()
    database.commit()

    return True


def update_homogenised_abundance(database, element, ion, cname,
    spectrum_filename_stub, abundance, uncertainty, num_lines, num_measurements,
    **kwargs):
    """
    Update the homogenised abundance for a species (element and ion) as measured
    in a spectrum (CNAME and spectrum filename stub).

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

    raise NotImplementedError


def update_homogenised_line_abundance(database, element, ion, cname, 
    spectrum_filename_stub, wavelength, abundance, uncertainty,
    num_measurements, node_bitmask=0):
    """
    Update the homogenised line abundance for a given CNAME, species (element
    and ion), and wavelength.

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

    return data.retrieve(database, """INSERT INTO homogenised_line_abundances
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


def match_node_line_abundances(database, element, ion, wavelength,
    column="scaled_abundance", tol=0.1, ignore_gaps=True):


    # Data.
    measurements = data.retrieve_table(database,
        """SELECT * FROM line_abundances WHERE element = %s AND ion = %s AND
        flags = 0 AND wavelength >= %s AND wavelength <= %s""",
        (element, ion, wavelength - tol, wavelength + tol))
    measurements = measurements.group_by(["spectrum_filename_stub"])

    nodes = sorted(set(measurements["node"]))
    
    # Build the matrix.
    X = np.nan * np.ones((len(nodes), len(measurements.groups)))
    for j, group in enumerate(measurements.groups):
        for row in group:
            X[nodes.index(row["node"]), j] = row[column]

    row_valid = np.all if ignore_gaps else np.any
    X = X[:, row_valid(np.isfinite(X), axis=0)]    
    return (np.ma.array(X, mask=~np.isfinite(X)), nodes)


def _homogenise_spectrum_line_abundances(database, element, ion, measurements,
    rho, nodes, column="scaled_abundance", e_column="e_abundance", **kwargs):

    # Note that by default we use scaled_abundance, not abundance, such that we
    # can apply any offsets or scaling as necessary *before* the homogenisation
    # begins.


    N_measurements = np.isfinite(measurements[column]).sum()

    if N_measurements == 0:
        return (np.nan, np.nan, 0, -1)

    # Problems to address later:
    assert np.all(measurements["upper_abundance"] == 0)

    abundance = np.nan * np.ones(len(nodes))
    u_abundance = np.nan * np.ones(len(nodes))
    for i, node in enumerate(nodes):
        if node not in measurements["node"]: continue
        match = (measurements["node"] == node)
        abundance[i] = measurements[column][match][0]
        u_abundance[i]  = measurements[e_column][match][0]

    mask = np.isfinite(abundance)
    abundance = abundance[mask]
    u_abundance = u_abundance[mask]

    # Fill in missing abundance uncertainties with the mean of the others.
    bad_values = ~np.isfinite(u_abundance) + (0 >= u_abundance)
    u_abundance[bad_values] = np.nanmean(u_abundance[~bad_values])
    if not np.any(np.isfinite(u_abundance)):
        logger.warn("Setting 0.2 dex uncertainty for {0} {1} line at {2}".format(
            element, ion, measurements["wavelength"][0]))
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

    assert np.isfinite(abundance_mean * abundance_sigma)

    # Place this information in the homogenised abundance table.
    unique_id = update_homogenised_line_abundance(database, element, ion,
        measurements["cname"][0], measurements["spectrum_filename_stub"][0],
        measurements["wavelength"][0], abundance_mean, abundance_sigma, N)
    
    return (abundance_mean, abundance_sigma, N, unique_id)


def homogenise_line_abundances(database, element, ion, cname, wavelength,
    rho, nodes, wavelength_tolerance=0.1, **kwargs):
    """
    Homogenise the line abundances for a given element, ion, star and wavelength
    (within some tolerance).

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

    :param wavelength_tolerance: [optional]
        The acceptable tolerance in wavelength for this element and ion.

    :type wavelength_tolerance:
        float
    """

    # Get all valid line abundances for this element, ion, cname.
    wavelength_tolerance = abs(wavelength_tolerance)
    if wavelength_tolerance > 0:   
        measurements = data.retrieve_table(database, """SELECT * 
            FROM line_abundances
            WHERE element = %s AND ion = %s AND cname = %s AND flags = 0
            AND wavelength >= %s AND wavelength <= %s ORDER BY node ASC""",
            (element, ion, cname, wavelength - wavelength_tolerance,
                wavelength + wavelength_tolerance))
    else:
        measurements = data.retrieve_table(database, """SELECT *
            FROM line_abundances WHERE element = %s AND ion = %s AND cname = %s
            AND flags = 0 AND wavelength = %s ORDER BY node ASC""",
            (element, ion, cname, wavelength))

    assert measurements is not None

    # Account for repeat spectra for a single star.
    results = {}
    measurements = measurements.group_by(["spectrum_filename_stub"])
    for i, group in enumerate(measurements.groups):
        stub = group["spectrum_filename_stub"][0]
        results[stub] = _homogenise_spectrum_line_abundances(database, element,
            ion, group, rho, nodes, **kwargs)

    return results


def homogenise_average_abundances(database, element, ion, cname):
    """
    Average the homogenised line abundances for a given star.
    """

    # Get the data for this star/spectrum filename stub.
    line_abundances = data.retrieve_table(database,
        """SELECT * FROM homogenised_line_abundances WHERE element = %s
        AND ion = %s AND cname = %s""", (element, ion, cname))
    assert line_abundances is not None

    line_abundances = line_abundances.group_by(["spectrum_filename_stub"])
    for i, group in enumerate(line_abundances.groups):

        # Calculate a weighted mean and variance.
        # TODO: currently ignoring the covariance between lines.
        inv_variance = 1.0/(group["e_abundance"]**2)
        weights = inv_variance/np.sum(inv_variance)

        mu = np.sum(group["abundance"] * weights)
        sigma = np.sqrt(np.sum(weights**2 * group["e_abundance"]**2))

        # Update the line abundance.
        N_lines = len(group)
        N_measurements = sum(line_abundances["num_measurements"])

        result = update_homogenised_abundance(database, element, ion,
            group["cname"][0], group["spectrum_filename_stub"][0],
            mu, sigma, N_lines, N_measurements)

        raise a
    raise a
    # Calculate a mean and std.dev.

    # Update this value into the homogenised_node_results table.


    raise NotImplementedError


def homogenise_species(database, element, ion, **kwargs):
    """
    Homogenise the line abundances for all stars for a given element and ion.

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
    """

    # Get the unique wavelengths.
    wavelengths = data.retrieve_column(database,
        """SELECT DISTINCT ON (wavelength) wavelength FROM line_abundances WHERE
        element = %s AND ion = %s AND flags = 0 ORDER BY wavelength ASC""",
        (element, ion), asarray=True)
    decimals = kwargs.pop("round_wavelengths", 1)
    if decimals >= 0: np.round(wavelengths, decimals, out=wavelengths)

    # Get the unique CNAMEs.
    cnames = data.retrieve_column(database,
        """SELECT DISTINCT ON (cname) cname FROM line_abundances
        WHERE element = %s AND ion = %s ORDER BY cname ASC""", (element, ion))

    # For each wavelength, approximate the covariance matrix then homogenise
    # this wavelength for all cnames.
    t_init = time()
    for i, wavelength in enumerate(set(wavelengths)):
        # Calculate the correlation coefficients for this line.
        X, nodes = match_node_line_abundances(database, element, ion, wavelength,
            **kwargs)
        rho = np.corrcoef(X)

        # For each cname, homogenise this line.
        for j, cname in enumerate(cnames):
            homogenise_line_abundances(database, element, ion, cname,
                wavelength, rho, nodes)

    # Need to commit before we can do the averaged results per star.
    database.commit()

    for j, cname in enumerate(cnames):
        homogenise_average_abundances(database, element, ion, cname)


    # Now produced homogenised results from all the available lines for this
    # species.

    # For each cname, just calculate a weighted average.
    # Include the line-to-line covariances?


    raise a

    # Get the matched line abundances for this species so we can estimate the
    # covariance matrix.

    # species
    # per cnames
    # per wavelength

    # Instead we should do:
    # per wavelength (get the covariance matrix)
    # then per cname
    


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

    import flags

    drop_table(database) # Drop the table if it already exists.
    create_table(database)

    """
    flags.create_line_abundance_flag(database, "line abundance is unphysically high")
    flags.create_line_abundance_flag(database, "i haven't had coffee yet")
    flags.create_line_abundance_flag(database, "my wallet is missing")
    flags.create_line_abundance_flag(database, "i love coffee yet")


    moo = flags.search_line_abundance_flags(database, "coffee")
    for k in moo:
        assert flags.flag_exists(database, int(k))


    flags.update_line_abundance_flag(database, [1, 2, 3], "id IN %s", (1, 2, 3, 4))
    
    """

    # It should always be that we do the following:

    # [ ] Mark line abundances that are obviously incorrect and should be removed
    # [ ] Calculate and apply any line/node-specific biases to the distributions
    # [X] The line abundances should be able to be homogenised per star/cname in
    #     this file.
    #

    # That means we need to be able to:
    # [X] Select line abundances based on any constraints.
    # [ ] Mark line abundances with some flag, with any constraints.
    # [X] Give multiple flags for lines.

    # [ ] Calculate and apply any line/node-specific biases to the distributions
    #     (this should only use abundances that are not marked)

    # The homogenised thing should just ignore anything with a mark > 0



    homogenise_species(database, "Si", 2)


    raise a

