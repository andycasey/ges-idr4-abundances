
""" Produce ensemble Mo 1 abundances from the Gaia-ESO Survey iDR4 WG11 data """

__author__ = 'Andy Casey <arc@ast.cam.ac.uk>'


import release

database, element, ion = ("arc", "Mo", 1)
ges = release.DataRelease(database)

# ULB report unphysically high values for Mo 1
flag_id = ges.flags.retrieve_or_create("Abundance is non-physical")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = '{1}' AND
        (abundance > 4)""".format(element, ion))

# ULB and MyGIsFOS are the only nodes that report Mo 1 abundances, and they do
# not use any common lines. The abundances show a large scatter below 3750 K and
# greater than ~5500 K.
flag_id = ges.flags.retrieve_or_create(
    "Abundance values show large scatter in this temperature range")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT cname, teff
        FROM node_results WHERE node like 'Lumba%') n
    ON (l.cname = n.cname AND element = '{0}' AND ion = '{1}'
    AND node LIKE 'ULB%' AND (teff < 3750 OR teff > 5500))""".format(
        element, ion))

# There are no differential biases for this element. The solar biases don't
# apply either because Mo 1 was not measured in the Sun by MyGIsFOS, and the
# solar temperature is outside of the valid range that we just applied.

# Calculate biases and apply them.
species_biases = ges.biases.solar(element, ion)
for node in species_biases:
    for wavelength, (bias, sigma, N) in species_biases[node].items():
        rows = ges.biases.apply_offset(element, ion, node, wavelength, -bias)

# Perform the homogenisation.
ges.homogenise.species(element, ion)

