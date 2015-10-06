
""" Produce ensemble Cu 1 abundances from the Gaia-ESO Survey iDR4 WG11 data """

__author__ = 'Andy Casey <arc@ast.cam.ac.uk>'


import release

database, element, ion = ("arc", "Cu", 1)
ges = release.DataRelease(database)


# CAUP showed large scatter for 5218.2 line
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1} AND
    node like 'CAUP%' AND wavelength < 5219""".format(element, ion))

# EPINARBO shows poor results for most EMP lines 
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this line for EMP stars")
num_rows = ges.flags.update([flag_id],
    """SELECT l.id FROM line_abundances l JOIN 
    (SELECT DISTINCT ON (cname) cname, feh from node_results) n 
    ON (l.cname = n.cname and l.element = '{0}' AND l.ion = {1}
        AND l.node LIKE 'EPINARBO%' AND n.feh < -2.5 and l.wavelength < 5700.3)"""\
    .format(element, ion))

# ULB bad for cool stars for certain lines.
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this star for cool stars")
num_rows = ges.flags.update([flag_id],
    """SELECT l.id FROM line_abundances l JOIN 
    (SELECT DISTINCT ON (cname) cname, teff from node_results) n 
    ON (l.cname = n.cname and l.element = '{0}' AND l.ion = {1}
        AND l.node LIKE 'ULB%' AND n.teff < 4000 
        AND ((l.wavelength > 5700 AND l.wavelength < 5701) OR 
        (l.wavelength < 5106)))"""\
    .format(element, ion))

# ULB bad for hot stars for one line
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this star for hot stars")
num_rows = ges.flags.update([flag_id],
    """SELECT l.id FROM line_abundances l JOIN 
    (SELECT DISTINCT ON (cname) cname, teff from node_results) n 
    ON (l.cname = n.cname and l.element = '{0}' AND l.ion = {1}
        AND l.node LIKE 'ULB%' AND n.teff > 6000 
        AND l.wavelength > 5220 AND l.wavelength < 5221)"""\
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
