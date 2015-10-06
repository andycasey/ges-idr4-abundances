
""" Produce ensemble Ca 1 abundances from the Gaia-ESO Survey iDR4 WG11 data """

__author__ = 'Andy Casey <arc@ast.cam.ac.uk>'


import release

database, element, ion = ("arc", "Ca", 1)
ges = release.DataRelease(database)

# LUMBA gives unphysical results.
# anything outside of 3.5, 8
flag_id = ges.flags.retrieve_or_create("Abundance is non-physical")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1}
    AND (abundance < 3.5 OR abundance > 8)""".format(element, ion))

# CAUP shows large scatter for all lines blueward of 5876 (mostly MR)
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1} AND
    node LIKE 'CAUP%' AND wavelength < 5876""".format(element, ion))

# CAUP shows large scatter for all lines redward of 6166 (mostly MR)
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1} AND
    node LIKE 'CAUP%' AND wavelength > 6166""".format(element, ion))

# EPINARBO shows scatter for EMP stars for 5582.0 line
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this line for EMP stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        feh FROM node_results) n ON (l.cname = n.cname AND l.element = '{0}'
        AND l.ion = {1} AND n.feh < -2.5 AND l.node LIKE 'EPINARBO%' AND
        l.wavelength > 5581.5 AND l.wavelength < 5582.5)""".format(element, ion))

# EPINARBO shows large scatter for 6122.2 line and 6163.8 line and 6169.0 line
# and for all lines redder than 6471.7 line (incl) until 6700.
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this line for EMP stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}'
        AND ion = {1} AND node LIKE 'EPINARBO%' AND
        (
            (wavelength > 6122 AND wavelength < 6123) OR 
            (wavelength > 6163 AND wavelength < 6164) OR
            (wavelength > 6471.0 AND wavelength < 6700)
        )""".format(element, ion))

# LUMBA shows large scatter for EMP stars for 5260.4 line
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this line for EMP stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        feh FROM node_results) n ON (l.cname = n.cname and l.element = '{0}'
        AND l.ion = {1} AND n.feh < -2.5 AND l.wavelength > 5260 AND
        l.wavelength < 5261 AND l.node LIKE 'Lumba%')""".format(element, ion))

# Lumba shows large scatter for MR stars at 5261.7
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        feh FROM node_results) n ON (l.cname = n.cname and l.element = '{0}'
        AND l.ion = {1} AND n.feh > 0 AND l.wavelength > 5261 AND
        l.wavelength < 5262 AND l.node LIKE 'Lumba%')""".format(element, ion))

# Lumba shows large scatter for EMP stars for 5867.6 line and 6156.0 lineABD 6166.4 line and 6455.6 line
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this line for EMP stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        feh FROM node_results) n ON (l.cname = n.cname and l.element = '{0}'
        AND l.ion = {1} AND n.feh < -2.5 AND l.node LIKE 'Lumba%'
        AND (
            (l.wavelength > 5867 AND l.wavelength < 5868) OR
            (l.wavelength > 6155.0 AND l.wavelength < 6156.5) OR
            (l.wavelength > 6166.0 AND l.wavelength < 6167) OR
            (l.wavelength > 6455.0 AND l.wavelength < 6456.0)
        )""".format(element, ion))

# Lumba shows large scatter for hot stars for 5601.3 and 5857.4 line and 
# 6102.7 and 6122.2 line and 6162.2 line
# Lumba is probably quite bad for all hot stars (>6250)
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this star for hot stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        teff FROM node_resutls) n ON (l.cname = n.cname AND l.element = '{0}'
        AND l.ion = {1} and n.teff > 6250 AND l.node LIKE 'Lumba%')""".format(
        element, ion))


# ULB is bad for EMP stars for 5260.4 line
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this star for EMP stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        feh FROM node_results) n ON (l.cname = n.cname AND l.element = '{0}'
        AND l.ion = {1} AND l.node LIKE 'ULB%' AND n.feh < -2.5
        AND l.wavelength > 5260 AND l.wavelength < 5261)""".format(
        element, ion))

# ULB shows large scatter for 5261.7 line
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1} AND
    node LIKE 'ULB%' AND wavelength > 5261 AND wavelength < 5262""".format(element, ion))

# ULB bad for cool stars in 5513.0 line
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this star for cool stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        teff FROM node_results) n ON (l.cname = n.cname AND l.element = '{0}'
        AND l.ion = {1} AND l.node LIKE 'ULB%' AND n.teff < 4000
        AND l.wavelength > 5512.5 AND l.wavelength < 5513.5)""".format(
        element, ion))

# ULB bad for 5582.0 line (large scatter) for cool (teff < 4800) MR (> -1) stars
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this star for cool stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        teff, feh FROM node_results) n ON (l.cname = n.cname AND l.element = '{0}'
        AND l.ion = {1} AND l.node LIKE 'ULB%' AND n.teff < 4000 AND n.feh > -1
        AND (
            (l.wavelength > 5581.5 AND l.wavelength < 5582.5)
        )""".format(element, ion))

# ULB bad for cool stars (<4000) and hot stars (>6000) in 5588.8 line abd 5857.4 line
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this star for cool stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        teff FROM node_results) n ON (l.cname = n.cname AND l.element = '{0}'
        AND l.ion = {1} AND l.node LIKE 'ULB%' AND n.teff < 4000
        AND (
            (l.wavelength > 5588 AND l.wavelength < 5589) OR 
            (l.wavelength > 5857 AND l.wavelength < 5858)
        )""".format(element, ion))

flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this star for hot stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        teff FROM node_results) n ON (l.cname = n.cname AND l.element = '{0}'
        AND l.ion = {1} AND l.node LIKE 'ULB%' AND n.teff > 6000
        AND (
            (l.wavelength > 5588 AND l.wavelength < 5589) OR 
            (l.wavelength > 5857 AND l.wavelength < 5858)
        )""".format(element, ion))

# ULB bad for metal-oor stars in 5867.6 line
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this star for EMP stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        feh FROM node_results) n ON (l.cname = n.cname AND l.element = '{0}'
        AND l.ion = {1} AND l.node LIKE 'ULB%' AND n.feh < -2.5
        AND l.wavelength > 5867 AND l.wavelength < 5868)""".format(
        element, ion))

# ULB bad for cool stars in 6102.7 line and 6122.2 line
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this star for cool stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        teff FROM node_results) n ON (l.cname = n.cname AND l.element = '{0}'
        AND l.ion = {1} AND l.node LIKE 'ULB%' AND n.teff < 4000
        AND l.wavelength > 6102.5 AND l.wavelength < 6103)""".format(
        element, ion))

# ULB bad for hot stars (>6000) in 6102.7 line and 6122.2 line
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this star for hot stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        teff FROM node_results) n ON (l.cname = n.cname AND l.element = '{0}'
        AND l.ion = {1} AND l.node LIKE 'ULB%' AND n.teff > 6000 AND
        (
            (l.wavelength > 6102 AND l.wavelength < 6103) OR
            (l.wavelength > 6122 AND l.wavelength < 6123)
        )""".format(element, ion))

# ULB bad for EMP stars in 6156.0 line
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this star for EMP stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        feh FROM node_results) n ON (l.cname = n.cname AND l.element = '{0}'
        AND l.ion = {1} AND l.node LIKE 'ULB%' AND n.feh < -2.5
        AND l.wavelength > 6155.5 AND l.wavelength < 6156.5)""".format(
        element, ion))

# ULB bad for cool stars for lines redder than 6156.0 
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this star for cool stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        teff FROM node_results) n ON (l.cname = n.cname AND l.element = '{0}'
        AND l.ion = {1} AND l.node LIKE 'ULB%' AND n.teff < 4000
        AND l.wavelength > 6156)""".format(
        element, ion))

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
