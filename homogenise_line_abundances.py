#!/usr/bin/python

""" Homogenise line abundances. """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"


import data
import numpy as np


def create_table(database):
    """
    Create a table for the homogenised line abundances.

    :param database:
        A PostgreSQL database connection.

    :type database:
        :class:`~psycopg2.connection`
    """

    # Superfluous columns removed:
    #object char(21) not null,
    #spectrum_filename_stub char(140) not null,

    cursor = database.cursor()
    cursor.execute("""CREATE TABLE homogenised_line_abundances (
        cname char(16) not null,
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
    cursor.close()
    database.commit()

    return True


def update_homogenised_line_abundance(database, element, ion, cname, wavelength,
    abundance, uncertainty, num_measurements, node_bitmask=0):
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

    cursor = database.cursor()
    cursor.execute("""INSERT INTO homogenised_line_abundances(cname, wavelength,
        element, ion, abundance, e_abundance, upper_abundance, num_measurements,
        node_bitmask) VALUES (%(cname)s, %(wavelength)s, %(element)s, %(ion)s,
        %(abundance)s, %(e_abundance)s, %(upper_abundance)s,
        %(num_measurements)s, %(node_bitmask)s) RETURNING id;""", {
            "cname": cname,
            "wavelength": wavelength,
            "element": element,
            "ion": ion,
            "abundance": abundance,
            "e_abundance": uncertainty,
            "upper_abundance": 0, # TODO
            "num_measurements": num_measurements,
            "node_bitmask": node_bitmask
    })
    unique_id = int(cursor.fetchone()[0])
    cursor.close()
    database.commit()

    return unique_id


def homogenise_line_abundances(database, element, ion, cname, wavelength,
    wavelength_tolerance=0.1):
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
            AND wavelength >= %s AND wavelength <= %s""", (element, ion, cname,
                wavelength - wavelength_tolerance,
                wavelength + wavelength_tolerance))
    else:
        measurements = data.retrieve_table(database, """SELECT *
            FROM line_abundances WHERE element = %s AND ion = %s AND cname = %s
            AND flags = 0 AND wavelength = %s""",
            (element, ion, cname, wavelength))

    assert measurements is not None

    """
    Now we have all line measurements for a given star, for a given wavelength
    (e.g., element and ion).
    """
    # Are there multiple spectrum_filename_stubs for this given CNAME?
    # (This occurs usually for benchmark stars, where there are multiple spectra
    #  for the same CNAME)
    spectra = set(measurements["spectrum_filename_stub"])

    # Note that we are using scaled_abundance, not abundance, such that we can
    # apply any offsets or scaling as necessary *before* the homogenisation
    # begins.

    N_spectra = len(spectra)
    N_measurements = np.isfinite(measurements["scaled_abundance"]).sum()

    # Problems to address later:
    try:

        assert np.all(measurements["upper_abundance"] == 0)
        assert N_spectra == 1
        assert N_measurements > 1

    except AssertionError:
        print("not doing anything because assert error")
        return (0, 0, 1, 1)

    """
    How do we homogenise the line abundances?

    TODO: For now we will just mean the abundances and take the std. dev.
    This means we are ignoring uncertainties provided by different nodes.
    TODO: We are also ignoring the covariances.
    # TODO: What if the abundance sigma is nan (1 measurement)

    """

    abundance_mean = np.nanmean(measurements["scaled_abundance"])
    abundance_sigma = np.nanstd(measurements["scaled_abundance"])

    # TODO: Create a node mask that tells us which nodes contributed here.

    # Place this information in the homogenised abundance table.
    result = update_homogenised_line_abundance(database, element, ion, cname,
        wavelength, abundance_mean, abundance_sigma, N_measurements)

    return (abundance_mean, abundance_sigma, N_measurements, result)


def homogenise_star_species(database, element, ion, cname, **kwargs):
    """
    Homogenise the line abundances for a given element, ion, and star (CNAME).

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
    """

    # Get the unique wavelengths.
    wavelengths = data.retrieve_column(database,
        """SELECT DISTINCT ON (wavelength) wavelength FROM line_abundances WHERE
        element = %s AND ion = %s AND cname = %s ORDER BY wavelength ASC""",
        (element, ion, cname), asarray=True)
    assert wavelengths is not None

    decimals = kwargs.pop("round_wavelengths", 1)
    if decimals >= 0: np.round(wavelengths, decimals, out=wavelengths)

    # Create a convenient lambda function to do the line abundances per lambda.
    h = lambda w: homogenise_line_abundances(database, element, ion, cname, w)

    return { w: h(w) for w in wavelengths }


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

    # Get the unique CNAMEs.
    cnames = data.retrieve_column(database,
        """SELECT DISTINCT ON (cname) cname FROM line_abundances
        WHERE element = %s AND ion = %s ORDER BY cname ASC""",
        (element, ion))
    assert cnames is not None
    from time import time

    t_init = time()
    results = []
    for cname in cnames[:1000]:
        results.append(homogenise_star_species(database, element, ion, cname, **kwargs))

    taken = (time() - t_init)
    print("Time taken for 10 stars in {0} {1}: {2:.0f} s".format(element, ion, taken))
    raise a




if __name__ == "__main__":

    import psycopg2 as pg

    database = pg.connect(dbname="arc")

    import flags

    drop_table(database) # Drop the table if it already exists.
    create_table(database)
    flags.create_line_abundance_flag(database, "line abundance is unphysically high")
    flags.create_line_abundance_flag(database, "i haven't had coffee yet")
    flags.create_line_abundance_flag(database, "my wallet is missing")
    flags.create_line_abundance_flag(database, "i love coffee yet")


    moo = flags.search_line_abundance_flags(database, "coffee")
    for k in moo:
        assert flags.flag_exists(database, int(k))


    flags.update_line_abundance_flag(database, [1, 2, 3], "id IN %s", (1, 2, 3, 4))
    

    raise a

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
