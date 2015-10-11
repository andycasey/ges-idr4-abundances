
""" Produce ensemble Ca 2 abundances from the Gaia-ESO Survey iDR4 WG11 data """

__author__ = 'Andy Casey <arc@ast.cam.ac.uk>'

import release

database, element, ion = ("arc", "V", 1)
ges = release.DataRelease(database)

# Note that V2 is only measured in a few stars by MYGISFOS, and the results look
# a little spurious, so we will not provide V2.

# Lots of nodes, lots of lines.

# Only MyGISFOS uses lines bluer than 4875
flag_id = ges.flags.retrieve_or_create("Not an ideal line for analysis")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE TRIM(element) = '{0}' AND ion = {1} AND
    wavelength < 4875""".format(element, ion))

# CAUP abundances are always systematically too high (compared to benchmarks)
# Suspect HFS/similar was not adopted.
# Looks like that CAUP removed their own results from the node results file.
flag_id = ges.flags.retrieve_or_create("Spurious results")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE TRIM(element) = '{0}' AND ion = {1} AND
    node LIKE 'CAUP%'""".format(element, ion))

# The 6216 line is no good.
flag_id = ges.flags.retrieve_or_create("Not an ideal line for analysis")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE TRIM(element) = '{0}' AND ion = {1} AND
    (wavelength > 6215 AND wavelength < 6216)""".format(element, ion))

# TODO: Lumba gives high results (>45)
flag_id = ges.flags.retrieve_or_create("Abundance is non-physical")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE TRIM(element) = '{0}' AND ion = {1} AND
    (abundance > 7 OR abundance < 0)""".format(element, ion))

# LUmba, MYGISFOS and Vilnius all agree quite well in mean abundances.
# EPINARBO seems to heavily underestimate V1 abundances compared tot hese nodes.

# 5698 line is no good.
flag_id = ges.flags.retrieve_or_create("Not an ideal line for analysis")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE TRIM(element) = '{0}' AND ion = {1} AND
    (wavelength > 5698 AND wavelength < 5699)""".format(element, ion))

# 5632 line is no good
flag_id = ges.flags.retrieve_or_create("Not an ideal line for analysis")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE TRIM(element) = '{0}' AND ion = {1} AND
    (wavelength > 5632 AND wavelength < 5633)""".format(element, ion))

# 4881 line is no good
flag_id = ges.flags.retrieve_or_create("Not an ideal line for analysis")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE TRIM(element) = '{0}' AND ion = {1} AND
    (wavelength > 4881 AND wavelength < 4882)""".format(element, ion))

# Calculate biases and apply them.
species_biases = ges.biases.differential(element, ion)
for node in species_biases:
    for wavelength, (bias, sigma, N) in species_biases[node].items():
        rows = ges.biases.apply_offset(element, ion, node, wavelength, -bias)

# Perform the homogenisation.
ges.homogenise.species(element, ion)

