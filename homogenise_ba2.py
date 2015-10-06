
""" Produce ensemble Ba 2 abundances from the Gaia-ESO Survey iDR4 WG11 data """

__author__ = 'Andy Casey <arc@ast.cam.ac.uk>'

import release

database, element, ion = ("arc", "Ba", 2)
ges = release.DataRelease(database)


# Remove EPINARBO measurements for feh < -2.5 for 4934.1 line
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this line for EMP stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        feh FROM node_results) n ON (l.cname = n.cname AND l.element = '{0}' AND
    l.ion = {1} AND l.wavelength > 4934 AND l.wavelength < 4935
    AND l.node LIKE 'EPINARBO%' AND n.feh < -2.5)""".format(element, ion))

# Remove LUMBA measurements for feh < -2.5 for 5853.7 line
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        feh FROM node_results) n ON (l.cname = n.cname AND l.element = '{0}' AND
    l.ion = {1} AND l.wavelength > 5853 AND l.wavelength < 5854
    AND l.node LIKE 'EPINARBO%' AND n.feh < -2.5)""".format(element, ion))

# Remove LUMBA measurements for teff < 4000 for 5853.7 line
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this star for cool stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        teff, feh FROM node_results) n ON (l.cname = n.cname AND
        l.element = '{0}' AND l.ion = {1} AND l.node LIKE 'Lumba%' AND
        l.wavelength > 5853 AND l.wavelength < 5854
        AND (n.teff < 4000 OR n.feh < -2.5))""".format(element, ion))


# EPINARBO is systematically higher than all other nodes, particularly for giants
# (Confirm that EPINARBO used HFS)

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