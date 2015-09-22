
''' Produce ensemble Si 1 abundances from the Gaia-ESO Survey iDR4 WG11 data '''

__author__ = 'Andy Casey <arc@ast.cam.ac.uk>'

import psycopg2 as pg
import biases, flags, homogenise, plot


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

# Use only a subset of trustworthy lines for all stars.
flag_id = flags.retrieve_or_create_line_abundance_flag(db,
    "Transition excluded by most nodes")
num_rows = flags.update_line_abundance_flag(db, [flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = '{1}' AND
    (wavelength < 4431.8 OR 
    (wavelength > 5772.1 AND wavelength < 6721.7))""".format(element, ion))

# Calculate biases and apply them.
species_biases = biases.differential_abundance_biases(db, element, ion)
for node in species_biases:
    for wavelength, (bias, sigma, N) in species_biases[node].items():
        rows = biases.apply_offset(db, element, ion, node, wavelength, -bias)

# Produce some figures.
fig = plot.differential_line_abundances(db, element, ion, scaled=True,
    ignore_flags=False)
fig.savefig("figures/{0}{1}/scaled-differential-line-abundances.png".format(
    element.upper(), ion))

fig = plot.compare_solar(db, element, ion, scaled=True, ignore_flags=False)
fig.savefig("figures/{0}{1}/scaled-compare-solar.png".format(element.upper(),
    ion))

raise a

# Perform the homogenisation.
result = homogenise.species(db, element, ion)



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
