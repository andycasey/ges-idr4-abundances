
""" Produce ensemble Cr 1 abundances from the Gaia-ESO Survey iDR4 WG11 data """

__author__ = 'Andy Casey <arc@ast.cam.ac.uk>'

import release

database, element, ion = ("arc", "Cr", 1)
ges = release.DataRelease(database)




# ULB shows very peculiar results.
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1} AND
    node like 'ULB%'""".format(element, ion))

raise NotImplementedError


# Calculate biases and apply them.
species_biases = ges.biases.differential(element, ion)
for node in species_biases:
    for wavelength, (bias, sigma, N) in species_biases[node].items():
        rows = ges.biases.apply_offset(element, ion, node, wavelength, -bias)

# Perform the homogenisation.
result = ges.homogenise.species(element, ion)

# Produce some figures.
fig = ges.plot.differential_line_abundances(element, ion, scaled=True)
fig.savefig("figures/{0}/scaled-differential-line-abundances.png".format(_))

fig = ges.plot.solar_comparison(element, ion, scaled=True)
fig.savefig("figures/{0}/scaled-compare-solar.png".format(_))

_ = "{0}{1}".format(element.upper(), ion)
fig = ges.plot.benchmark_line_abundances(element, ion, "benchmarks.yaml",
    scaled=True)
fig.savefig("figures/{0}/scaled-benchmark-line-abundances.png".format(_))
