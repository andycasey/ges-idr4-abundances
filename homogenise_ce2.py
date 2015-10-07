
""" Produce ensemble Ce 2 abundances from the Gaia-ESO Survey iDR4 WG11 data """

__author__ = 'Andy Casey <arc@ast.cam.ac.uk>'

import release

database, element, ion = ("arc", "Ce", 2)
ges = release.DataRelease(database)

# CAUP Ce 2 abundances for the 4773.9 line are mostly exactly zero. There must
# be a bug in their code.
flag_id = ges.flags.retrieve_or_create("Error suspected: abundance is frequently zero")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = '{1}' AND
        node LIKE 'CAUP%' AND wavelength < 4774""".format(element, ion))

# MyGIsFOS uses lines that none of the other nodes use.

# 4984 is no good.
# 4773.9 looks no good
# The lines 4984 and 4773.9 show peculiarities between nodes, and a large scatter.
flag_id = ges.flags.retrieve_or_create(
    "Line shows peculiarities between nodes, and a large overall scatter")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = '{1}'
    AND abundance != 'NaN' AND flags = 0 AND
    ((wavelength > 4773.8 AND wavelength < 4774) OR
    (wavelength > 4984.3 AND wavelength < 4984.5))""".format(
        element, ion))


# ULB measurements below teff < 3750
flag_id = ges.flags.retrieve_or_create(
    "Abundance values show large scatter in this temperature range")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT cname, teff
        FROM node_results WHERE node like 'Lumba%') n
    ON (l.cname = n.cname AND element = '{0}' AND ion = '{1}'
    AND node LIKE 'ULB%' AND (teff < 3750))""".format(
        element, ion))

# Calculate biases and apply them.
species_biases = ges.biases.differential(element, ion)
for node in species_biases:
    for wavelength, (bias, sigma, N) in species_biases[node].items():
        rows = ges.biases.apply_offset(element, ion, node, wavelength, -bias)

# Perform the homogenisation.
ges.homogenise.species(element, ion)

# Produce some figures.
fig = ges.plot.differential_line_abundances(element, ion, scaled=True)
fig.savefig("figures/{0}{1}/scaled-differential-line-abundances.png".format(
    element.upper(), ion))

fig = ges.plot.solar_comparison(element, ion, scaled=True)
fig.savefig("figures/{0}{1}/scaled-compare-solar.png".format(element.upper(),
    ion))
