
""" Produce ensemble Al 1 abundances from the Gaia-ESO Survey iDR4 WG11 data """

__author__ = 'Andy Casey <arc@ast.cam.ac.uk>'


import release

database, element, ion = ("arc", "Al", 1)
ges = release.DataRelease(database)


# Flag any lines that should be discarded, from specific nodes?
flag_id = ges.flags.retrieve_or_create("Abundance is non-physical")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = '{1}' AND
        (abundance > 8)""".format(element, ion))

# The Al 1 line at 6696.2 from MyGIsFOS is probably a duplicate of 6696.0
flag_id = ges.flags.retrieve_or_create("Line is probably a duplicate print of Al 1 at 6696.0")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = '{1}'
    AND node ILIKE 'MyGIsFOS%' AND wavelength > 6696.1 AND wavelength < 6696.2
    AND flags = 0""".format(element, ion))

# The Al 1 line at 6696.8 is not great for analysis.
flag_id = ges.flags.retrieve_or_create(
    "Not an ideal line for analysis")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = '{1}' AND 
    wavelength > 6696.7 AND wavelength < 6696.9 AND abundance != 'NaN'""".format(
        element, ion))

# ULB does not do well for 5557.1 for cool stars
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this star for cool stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname, 
        teff FROM node_results) n ON (n.cname = l.cname AND l.element = '{0}' 
    and l.ion = {1} and l.node LIKE 'ULB%' AND l.wavelength < 5557.2 AND n.teff < 4000)""".format(element, ion))

# Lumba does not do well for any metal-poor stars.
flag_id = ges.flags.retrieve_or_create("Large scatter seen in this line for EMP stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        feh FROM node_results) n ON (n.cname = l.cname AND l.element = '{0}'
    AND l.ion = {1} and l.node LIKE 'Lumba%' and n.feh < -2.5)""".format(element, ion))

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
