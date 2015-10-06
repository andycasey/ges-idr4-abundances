
""" Produce ensemble Si 1 abundances from the Gaia-ESO Survey iDR4 WG11 data """

__author__ = 'Andy Casey <arc@ast.cam.ac.uk>'

import release

database, element, ion = ("arc", "Si", 1)
ges = release.DataRelease(database)

# Flag any lines that should be discarded, from specific nodes?
flag_id = ges.flags.retrieve_or_create("Abundance is non-physical")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = '{1}' AND
        (abundance > 9 OR abundance < 5)""".format(element, ion))

# Use only a subset of trustworthy lines for all stars.
flag_id = ges.flags.retrieve_or_create("Transition excluded by most nodes")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = '{1}' AND
    (wavelength < 4431.8 OR 
    (wavelength > 5772.1 AND wavelength < 6721.7))""".format(element, ion))

# Calculate biases and apply them.
species_biases = ges.biases.differential(element, ion)
for node in species_biases:
    for wavelength, (bias, sigma, N) in species_biases[node].items():
        rows = ges.biases.apply_offset(element, ion, node, wavelength, -bias)

# Perform the homogenisation.
result = ges.homogenise.species(element, ion)

# Produce some figures.
fig = ges.plot.differential_line_abundances(element, ion, scaled=True)
fig.savefig("figures/{0}{1}/scaled-differential-line-abundances.png".format(
    element.upper(), ion))

fig = ges.plot.solar_comparison(element, ion, scaled=True)
fig.savefig("figures/{0}{1}/scaled-compare-solar.png".format(element.upper(),
    ion))
