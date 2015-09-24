
''' Produce ensemble Si 2 abundances from the Gaia-ESO Survey iDR4 WG11 data '''

__author__ = 'Andy Casey <arc@ast.cam.ac.uk>'


import release


element, ion = ("Si", 2)
ges = release.DataRelease("arc")


# Flag any lines that should be discarded, from specific nodes?
flag_id = ges.flags.retrieve_or_create("Abundance is non-physical")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = '{1}' AND
        (abundance > 9 OR abundance < 5)""".format(element, ion))
flag_id = ges.flags.retrieve_or_create(
    "Weak suspicious measurement based on comparisons to all other nodes")
num_rows2 = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = '{1}'
    AND node LIKE 'Lumba%' AND abundance < 5.8 AND flags = 0""".format(element, ion))

# Calculate biases and apply them.
species_biases = ges.biases.differential(element, ion)
for node in species_biases:
    for wavelength, (bias, sigma, N) in species_biases[node].items():
        rows = ges.biases.apply_offset(element, ion, node, wavelength, -bias)

# Produce some figures.
fig = ges.plot.differential_line_abundances(element, ion, scaled=True,
    ignore_flags=False)
fig.savefig("figures/{0}{1}/scaled-differential-line-abundances.png".format(
    element.upper(), ion))

fig = ges.plot.solar_comparison(element, ion, scaled=True, ignore_flags=False)
fig.savefig("figures/{0}{1}/scaled-compare-solar.png".format(element.upper(),
    ion))

# Perform the homogenisation.
result = ges.homogenise.species(element, ion)


