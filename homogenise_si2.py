
''' Produce ensemble Si 2 abundances from the Gaia-ESO Survey iDR4 WG11 data '''

__author__ = 'Andy Casey <arc@ast.cam.ac.uk>'

import psycopg2 as pg
import flags


db = pg.connect(dbname="arc")

element, ion = ("Si", 2)


# Flag any lines that should be discarded, from specific nodes?
#flag_id = flags.create_line_abundance_flag(db, "")
#flags.update_line_abundance_flag(db, [flag_id], )
#flags.update_line_abundance_flag(database, [flag_id], where, values=None):

# [EPINARBO] Any lines with log(Si 2) > 9 are unphysical (check they a)



# Calculate any biases for any lines, for any nodes?

# Apply biases.


# Perform the homogenisation.
homogenise_species(database, element, ion)

# Produce comparison figures.

# - node-to-node scatter plot
# - benchmark comparison
# - bensby comparison
# - 