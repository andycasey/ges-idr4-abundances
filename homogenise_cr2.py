
""" Produce ensemble Cr 2 abundances from the Gaia-ESO Survey iDR4 WG11 data """

__author__ = 'Andy Casey <arc@ast.cam.ac.uk>'

import release

database, element, ion = ("arc", "Cr", 2)
ges = release.DataRelease(database)

# 18 lines from 5 nodes

# 5305.8 and 5502.1 line very poor for Lumba in the sun.
# all nodes systematically high for 5420.9 line wrt to the sun

# remove cool (<4500) or EMP (<-2.5) stars for 5246.8 line for ULB
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this star for cool stars")

num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (select distinct on (cname) cname,
        teff, feh FROM node_results) n ON (l.cname = n.cname AND 
        l.element = '{0}' AND l.ion = {1} AND l.node LIKE 'ULB%' AND
        (n.feh < -2.5 OR n.teff < 4500) AND
        (l.wavelength > 5246 AND l.wavelength < 5247))""".format(
    element, ion))

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
