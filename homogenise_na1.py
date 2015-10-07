
""" Produce ensemble Na 1 abundances from the Gaia-ESO Survey iDR4 WG11 data """

__author__ = 'Andy Casey <arc@ast.cam.ac.uk>'


import release

database, element, ion = ("arc", "Na", 1)
ges = release.DataRelease(database)

# Non-physical errors
flag_id = ges.flags.retrieve_or_create("Abundance error is non-physical")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE TRIM(element) = '{0}' AND ion = '{1}' AND
        e_abundance > 1""".format(element, ion))

# Lines bluer than 4982.8 are generally only measured in a few stars by one node
flag_id = ges.flags.retrieve_or_create(
    "Not an ideal line for analysis")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = '{1}'
    AND wavelength < 4982""".format(element, ion))

# The 5890.0 and 5895.9 lines show a lot of scatter and irregularities.
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = '{1}'
    AND wavelength > 5789 AND wavelength < 5896""".format(element, ion))

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
