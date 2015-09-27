
""" Produce ensemble La 2 abundances from the Gaia-ESO Survey iDR4 WG11 data """

__author__ = 'Andy Casey <arc@ast.cam.ac.uk>'

import release

database, element, ion = ("arc", "La", 2)
ges = release.DataRelease(database)

# Up to 13 lines mesaured from 3 nodes (ULB, EPINARBO, and MYGIsFOS)
# ULB gives odd results: clip below -3 (and perhaps above +3)
flag_id = ges.flags.retrieve_or_create("Abundance is non-physical")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = '{1}' AND
        (abundance < -3)""".format(element, ion))

# Only the 5303.5 and 6390.5 lines are common to more than one node
# (EPINARBO AND ULB) and they both disagree with each other by ~0.3 dex.

flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this line for EMP stars")

# EPINARBO @ 5303 FOR FEH < -2
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT (cname) cname, feh
        FROM node_results) n ON (
        l.element = '{0}' AND l.ion = '{1}' AND l.node LIKE 'EPINARBO%'
        AND l.cname = n.cname and n.feh < -2 AND
        (l.wavelength > 5302 AND l.wavelength < 5304))""".format(element, ion))

# EPINARBO @ 4804 for FEH < -2.5
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT (cname) cname, feh
        FROM node_results) n ON (
        l.element = '{0}' AND l.ion = '{1}' AND l.node LIKE 'EPINARBO%'
        AND l.cname = n.cname and n.feh < -2.5 AND
        (l.wavelength > 4803 AND l.wavelength < 4805))""".format(element, ion))

# EPINARBO @ 5936.2 for FEH < -2
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT (cname) cname, feh
        FROM node_results) n ON (
        l.element = '{0}' AND l.ion = '{1}' AND l.node LIKE 'EPINARBO%'
        AND l.cname = n.cname and n.feh < -2 AND
        (l.wavelength > 5936 AND l.wavelength < 5937))""".format(element, ion))

# 4921.8 line looks no good
flag_id = ges.flags.retrieve_or_create("Not an ideal line for analysis")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = '{1}' AND
        (wavelength > 4921 AND wavelength < 4922)""".format(element, ion))

# 5123 line looks no good
flag_id = ges.flags.retrieve_or_create("Not an ideal line for analysis")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = '{1}' AND
        (wavelength > 5122.5 AND wavelength < 5123.5)""".format(element, ion))


# Note there is a strong positive correlation between LA 2 abundance and logg

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
