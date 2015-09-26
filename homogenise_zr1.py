
""" Produce ensemble Zr 1 abundances from the Gaia-ESO Survey iDR4 WG11 data """

__author__ = 'Andy Casey <arc@ast.cam.ac.uk>'


import release

database, element, ion = ("arc", "Zr", 1)
ges = release.DataRelease(database)


# Lumba and ULB give some log(X) abundances greater than 30.
flag_id = ges.flags.retrieve_or_create("Abundance is non-physical")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = '{1}' AND
        (abundance > 6 OR abundance < 0)""".format(element, ion))

# ULB Shows unusual results for the 4784.9 and 4828.0 line
flag_id = ges.flags.retrieve_or_create(
    "Line showed peculiar abundances with respect to stellar parameters")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = '{1}'
    AND node ILIKE 'ULB%' AND wavelength < 4785 AND abundance != 'NaN'""".format(
        element, ion))

# EPINARBO is the only node to use 4805.9, 4815.0, 5046.6, and 6445.7 lines
# The abundances seem to look OK.

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
