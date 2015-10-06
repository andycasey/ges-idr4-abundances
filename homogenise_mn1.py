
""" Produce ensemble Mn 1 abundances from the Gaia-ESO Survey iDR4 WG11 data """

__author__ = 'Andy Casey <arc@ast.cam.ac.uk>'


import release

database, element, ion = ("arc", "Mn", 1)
ges = release.DataRelease(database)

# Lumba routinely gives abundances of like 30+
flag_id = ges.flags.retrieve_or_create("Abundance is non-physical")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1}
    AND (abundance < 0 OR abundance > 8)""".format(element, ion))

# Lumba is bad for EMP stars for all lines except 5457.5 line
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this line for EMP stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        feh FROM node_results) n ON (l.element = '{0}' AND l.ion = {1} AND
        l.cname = n.cname AND l.node LIKE 'Lumba%' AND n.feh < -2.5 AND
        NOT (l.wavelength > 5457 AND l.wavelength < 5458))""".format(element, ion))

# Lumba shows a lot of scatter for the 5394.6 line
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1} AND
    node LIKE 'Lumba%' AND wavelength > 5394 AND wavelength < 5395""".format(
        element, ion))

# Lumba is bad for cool stars for 5117.9 line
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this star for cool stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        teff FROM node_results) n ON (l.element = '{0}' AND l.ion = {1} AND
        l.cname = n.cname AND l.node LIKE 'Lumba%' AND n.teff < 4000 AND
        l.wavelength > 5117 AND l.wavelength < 5118)""".format(element, ion))

# ULB for 4783.4 is wrong (large scatter)
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1} AND
    node LIKE 'ULB%' AND wavelength > 4783 AND wavelength < 4784""".format(
        element, ion))

# ULB is bad for EMP lines for 5004.9 line 
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this line for EMP stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        feh FROM node_results) n ON (l.cname = n.cname AND l.element = '{0}'
        AND l.ion = {1} AND n.feh < -2.5 AND l.node LIKE 'ULB%'
        AND l.wavelength > 5004 AND l.wavelength < 5005)""".format(element, ion))

# ULB is bad for cool stars for 5117.9 line
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this star for cool stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        teff FROM node_results) n ON (l.element = '{0}' AND l.ion = {1} AND
        l.cname = n.cname AND l.node LIKE 'ULB%' AND n.teff < 4000 AND
        l.wavelength > 5117 AND l.wavelength < 5118)""".format(element, ion))

# ULB is bad for EMP stars for 5255.4 line
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this line for EMP stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        feh FROM node_results) n ON (l.cname = n.cname AND l.element = '{0}'
        AND l.ion = {1} AND n.feh < -2.5 AND l.node LIKE 'ULB%'
        AND l.wavelength > 5255 AND l.wavelength < 5256)""".format(element, ion))

# ULB shows spurious results for 5394.6 line for metal-rich stars
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        feh FROM node_results) n ON (l.element = '{0}' AND l.ion = {1} AND
        l.cname = n.cname AND l.node LIKE 'ULB%' AND n.feh > 0 AND
        l.wavelength > 5394 AND l.wavelength < 5395)""".format(element, ion))

# ULB is bad for hot (>6000) stars at 5407.5 line
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this star for hot stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        teff FROM node_results) n ON (l.element = '{0}' AND l.ion = {1} AND
        l.cname = n.cname AND l.node LIKE 'ULB%' AND n.teff > 6000 AND
        l.wavelength > 5407 AND l.wavelength < 5408)""".format(element, ion))

# ULB is bad for 5432.5 line for cool stars (<4800)
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this star for cool stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        teff FROM node_results) n ON (l.element = '{0}' AND l.ion = {1} AND
        l.cname = n.cname AND l.node LIKE 'ULB%' AND n.teff < 4800 AND
        l.wavelength > 5432 AND l.wavelength < 5433)""".format(element, ion))

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
