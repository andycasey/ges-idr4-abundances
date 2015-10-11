
""" Produce ensemble Sr 1 abundances from the Gaia-ESO Survey iDR4 WG11 data """

__author__ = 'Andy Casey <arc@ast.cam.ac.uk>'


import release

database, element, ion = ("arc", "Sr", 1)
ges = release.DataRelease(database)

# ULB is bad for all hot (>5250) stars
flag_id = ges.flags.retrieve_or_create("Large scatter seen in this star for hot stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        teff FROM node_results) n ON (l.cname = n.cname AND l.element = '{0}'
        AND l.ion = {1} AND l.node LIKE 'ULB%' AND n.teff > 5250)""".format(
        element, ion))

# Calculate biases and apply them.
species_biases = ges.biases.solar(element, ion)
for node in species_biases:
    for wavelength, (bias, sigma, N) in species_biases[node].items():
        rows = ges.biases.apply_offset(element, ion, node, wavelength, -bias)

# Perform the homogenisation.
ges.homogenise.species(element, ion)

