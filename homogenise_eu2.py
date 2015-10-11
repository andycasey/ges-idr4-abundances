
""" Produce ensemble Eu 2 abundances from the Gaia-ESO Survey iDR4 WG11 data """

__author__ = 'Andy Casey <arc@ast.cam.ac.uk>'

import release

database, element, ion = ("arc", "Eu", 2)
ges = release.DataRelease(database)

# 5 lines from 3 nodes (EPINARBO, Lumba, ULB)

# Abundances greater than 4 should be clipped.
# The 5818 line is only used by EPINARBO, and only very occasionally.

# EPINARBO does not supply 6645 line for most stars below [Fe/H] < -1.5, which
# is why some of the plots look like there is a 'clump' below the normal Eu 2
# distribution

# ULB and Lumba abundances for Solar look reasonable. No EPINARBO Solar results
# available.

# EPINARBO results look systematically higher (~0.2 dex) than the other nodes.

# The 6049 and 6173 lines also show a lot of scatter.
flag_id = ges.flags.retrieve_or_create("Abundance is non-physical")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = '{1}' AND
        (abundance > 4)""".format(element, ion))

flag_id = ges.flags.retrieve_or_create(
    "Not an ideal line for analysis")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = '{1}'
    AND wavelength < 6437""".format(element, ion))


# Solar biases available for most nodes. Differential biases available for 1/2
# of the remaining lines.

# Calculate biases and apply them.
species_biases = ges.biases.differential(element, ion)
for node in species_biases:
    for wavelength, (bias, sigma, N) in species_biases[node].items():
        rows = ges.biases.apply_offset(element, ion, node, wavelength, -bias)

# Perform the homogenisation.
ges.homogenise.species(element, ion)

