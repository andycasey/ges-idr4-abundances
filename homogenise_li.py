
""" Produce ensemble Li 1 abundances from the Gaia-ESO Survey iDR4 WG11 data """

__author__ = 'Andy Casey <arc@ast.cam.ac.uk>'


import release

database, element, ion = ("arc", "Li", 1)
ges = release.DataRelease(database)

# Extract OACT results.


# Set all wavelengths to the same value.
ges.release.execute("""UPDATE line_abundances SET wavelength = 6707.8 WHERE
    element = '{0}' AND ion = {1}""".format(element, ion))
ges.release.commit()

# Need OACT results
# Need to deal with upper limits from EPINARBO gracefully
raise NotImplementedError


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
