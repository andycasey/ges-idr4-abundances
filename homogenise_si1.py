
''' Produce ensemble Si 1 abundances from the Gaia-ESO Survey iDR4 WG11 data '''

__author__ = 'Andy Casey <arc@ast.cam.ac.uk>'

import psycopg2 as pg
import biases, flags, homogenise


kwds = {
    "host": "/tmp/",
    "dbname": "arc",
    "user": "arc",
    "password": "password"
}

db = pg.connect(**kwds)

element, ion = ("Si", 1)

# Flag any lines that should be discarded, from specific nodes?
flag_id = flags.retrieve_or_create_line_abundance_flag(db,
    "Abundance is non-physical")
num_rows = flags.update_line_abundance_flag(db, [flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = '{1}' AND
        (abundance > 9 OR abundance < 5)""".format(element, ion))

# Any LUMBA Si 2 measurements below X (or REW before Y) should be flagged.
# They look 
flag_id = flags.retrieve_or_create_line_abundance_flag(db, 
    "Weak suspicious measurement based on comparisons to all other nodes")
num_rows2 = flags.update_line_abundance_flag(db, [flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = '{1}'
    AND node LIKE 'Lumba%' AND abundance < 6 AND flags = 0""".format(element, ion))


# Calculate biases and apply them.
solar_biases = biases.solar_abundance_biases(db, element, ion)
for node in solar_biases:
    for wavelength, (bias, sigma, N) in solar_biases[node].items():
        rows = biases.apply_offset(db, element, ion, node, wavelength, -bias)


# Produce some figures.
import plot

#fig = plot.differential_line_abundances(db, element, ion, scaled=True,
#    ignore_flags=True)

#raise a

# Perform the homogenisation.
from time import time
t = time()

result = homogenise.species(db, element, ion)
taken = time() - t
raise a
# Produce comparison figures.

# - node-to-node scatter plot
# - benchmark comparison
# - bensby comparison
# - 

raise a 
