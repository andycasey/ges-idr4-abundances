
""" Produce ensemble Ca 2 abundances from the Gaia-ESO Survey iDR4 WG11 data """

__author__ = 'Andy Casey <arc@ast.cam.ac.uk>'

import release

database, element, ion = ("arc", "Ca", 2)
ges = release.DataRelease(database)

# 3 lines from 3 nodes

# Lumba reports some unusual measurements.
flag_id = ges.flags.retrieve_or_create("Abundance is non-physical")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = '{1}' AND
        (abundance > 7.5 OR abundance < 3.5)""".format(element, ion))

# The 5339 line is always under-abundant, it seems.
flag_id = ges.flags.retrieve_or_create("Not an ideal line for analysis")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1} AND
    wavelength > 5339 AND wavelength < 5340""")

# Lumba is bad for metal-poor stars
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this line for EMP stars")
num_rows = ges.flags.update([flag_id],
    """SELECT l.id FROM line_abundances l JOIN 
    (SELECT DISTINCT ON (cname) cname, feh from node_results) n 
    ON (l.cname = n.cname and l.element = '{0}' AND l.ion = {1}
        AND l.node LIKE 'Lumba%' AND n.feh < -2.5)""".format(element, ion))

# ULB for 5001.5 line for metallicity > -0.25
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this line for ULB metal-rich stars")
num_rows = ges.flags.update([flag_id],
    """SELECT l.id FROM line_abundances l JOIN 
    (SELECT DISTINCT ON (cname) cname, feh from node_results) n 
    ON (l.cname = n.cname and l.element = '{0}' AND l.ion = {1}
        AND l.node LIKE 'ULB%' AND n.feh > -0.25 AND l.wavelength < 5002)"""\
    .format(element, ion))

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
