
""" Produce ensemble Pr 2 abundances from the Gaia-ESO Survey iDR4 WG11 data """

__author__ = 'Andy Casey <arc@ast.cam.ac.uk>'

import release

database, element, ion = ("arc", "Pr", 2)
ges = release.DataRelease(database)

# Only 4 lines present from two nodes (ULB and MyGIsFOS):
# 4222.9: MyGIsFOS only
# 4510.2: MyGIsFOS only
# 5259.7: ULB only 
# 5322.8: ULB only

# The 4510.2 line is not measured in the Sun.
# The ULB values for the Sun are reasonable.
# MyGIsFOS only measured Pr 2 in the Sun once, yielding ~ -0.4 dex off Solar.

# ULB gives odd results: clip below -3 (and perhaps above +3)
flag_id = ges.flags.retrieve_or_create("Abundance is non-physical")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = '{1}' AND
        (abundance < -3)""".format(element, ion))

# No differential biases available. Solar biases are available for most lines.
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
fig.savefig("figures/{0}{1}/scaled-differential-line-abundances.png".format(
    element.upper(), ion))

fig = ges.plot.solar_comparison(element, ion, scaled=True)
fig.savefig("figures/{0}{1}/scaled-compare-solar.png".format(element.upper(),
    ion))
