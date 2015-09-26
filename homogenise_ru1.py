
""" Produce ensemble Ru 1 abundances from the Gaia-ESO Survey iDR4 WG11 data """

__author__ = 'Andy Casey <arc@ast.cam.ac.uk>'

import release

database, element, ion = ("arc", "Ru", 1)
ges = release.DataRelease(database)

# We only have one line (4869.2) from one node (ULB).
# Just remove obviously spurious results.
flag_id = ges.flags.retrieve_or_create("Abundance is non-physical")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = '{1}' AND
        (abundance > 7 or abundance < -5)""".format(element, ion))

# No solar or differential biases available.
"""
# Calculate biases and apply them.
species_biases = ges.biases.differential(element, ion)
for node in species_biases:
    for wavelength, (bias, sigma, N) in species_biases[node].items():
        rows = ges.biases.apply_offset(element, ion, node, wavelength, -bias)
"""

# Perform the homogenisation.
ges.homogenise.species(element, ion)

# Produce some figures.
fig = ges.plot.differential_line_abundances(element, ion, scaled=True)
if fig is not None:
    fig.savefig("figures/{0}{1}/scaled-differential-line-abundances.png".format(
        element.upper(), ion))

fig = ges.plot.solar_comparison(element, ion, scaled=True)
if fig is not None:
    fig.savefig("figures/{0}{1}/scaled-compare-solar.png".format(element.upper(),
        ion))
