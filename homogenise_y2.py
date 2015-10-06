
""" Produce ensemble Y 2 abundances from the Gaia-ESO Survey iDR4 WG11 data """

__author__ = 'Andy Casey <arc@ast.cam.ac.uk>'

import release

database, element, ion = ("arc", "Y", 2)
ges = release.DataRelease(database)

# 21 lines from 6 nodes.

# Vilnius line abundances seem systematically lower than other nodes.
# MyGISfos line abundances are all systematically low for the Sun...

# 4854.9 line looks OK, but needs differential/solar offsets
# 4883.7 line looks OK (only EPINARBO)
# 4900.1 looks OK, modulo some offset






# remove 4204.7 line
# remove 4398.0 line
flag_id = ges.flags.retrieve_or_create("Not an ideal line for analysis")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE TRIM(element) = '{0}' AND ion = {1}
        AND (wavelength > 4204 AND wavelength < 4205)""".format(element, ion))

flag_id = ges.flags.retrieve_or_create("Not an ideal line for analysis")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE TRIM(element) = '{0}' AND ion = {1}
        AND (wavelength > 4397.5 AND wavelength < 4398.5)""".format(element, ion))

# remove 5119.1 and 5200.4 line for epinarbo for feh < -2.5
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this line for EMP stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (select distinct on (cname) cname,
        feh FROM node_results) n ON (l.cname = n.cname AND 
        TRIM(l.element) = '{0}' AND l.ion = {1} AND l.node LIKE 'EPINARBO%' AND
        n.feh < -2.5 AND (
            (l.wavelength > 5119 AND l.wavelength < 5120) OR
            (l.wavelength > 5199.5 AND l.wavelength < 5201.5)))""".format(
    element, ion))

# remove ulb and lumba for 4982.1 line for stars with feh < -2.5
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this line for EMP stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (select distinct on (cname) cname,
        feh FROM node_results) n ON (l.cname = n.cname AND 
        TRIM(l.element) = '{0}' AND l.ion = {1} AND 
        (l.node LIKE 'ULB%' OR l.node LIKE 'Lumba%')
        AND n.feh < -2.5 AND 
        (l.wavelength > 4982 AND l.wavelength < 4983))""".format(
    element, ion))

# remove ulb data for 5402.8 line because they do bad with cool stars
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this star for cool stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (select distinct on (cname) cname,
        feh FROM node_results) n ON (l.cname = n.cname AND 
        TRIM(l.element) = '{0}' AND l.ion = {1} AND 
        l.node LIKE 'ULB%' AND 
        (l.wavelength > 5402 AND l.wavelength < 5403))""".format(
    element, ion))

# remove ulb data for 5544.6 line for teff < 4000 K
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this star for cool stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (select distinct on (cname) cname,
        feh, teff FROM node_results) n ON (l.cname = n.cname AND 
        TRIM(l.element) = '{0}' AND l.ion = {1} AND 
        l.node LIKE 'ULB%' AND n.teff < 4000 AND 
        (l.wavelength > 5544 AND l.wavelength < 5545))""".format(
    element, ion))


# EPINARBO Bad for EMP stars at 5544.6 and 5546.0 line and 4982.1 line and 5119.1 line and 5123.2 line
# and 5200.4 line and 5289.8 line and 5402.8 line and 5473.4 line
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this line for EMP stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (select distinct on (cname) cname,
        feh FROM node_results) n ON (l.cname = n.cname AND 
        TRIM(l.element) = '{0}' AND l.ion = {1} AND l.node LIKE 'EPINARBO%'
        AND n.feh < -2.5 AND (
            (l.wavelength > 5544 AND l.wavelength < 5545) OR
            (l.wavelength > 5545.5 AND l.wavelength < 5546.5) OR
            (l.wavelength > 4982 AND l.wavelength < 4983) OR
            (l.wavelength > 5119 AND l.wavelength < 5120) OR
            (l.wavelength > 5123 AND l.wavelength < 5124) OR
            (l.wavelength > 5200 AND l.wavelength < 5201) OR
            (l.wavelength > 5289 AND l.wavelength < 5290) OR
            (l.wavelength > 5402 AND l.wavelength < 5403) OR
            (l.wavelength > 5473 AND l.wavelength < 5474)
        ))""".format(element, ion))

# Lumba bad for 4982.1 line for EMP lines and 5544.6 line
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this line for EMP stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (select distinct on (cname) cname,
        feh FROM node_results) n ON (l.cname = n.cname AND 
        TRIM(l.element) = '{0}' AND l.ion = {1} AND 
        l.node LIKE 'Lumba%'
        AND n.feh < -2.5 AND 
        (
            (l.wavelength > 5544 AND l.wavelength < 5545)
        ))""".format(
    element, ion))

# Lumba shows large scatter for 5728.9 line
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE TRIM(element) = '{0}' AND ion = {1}
    AND node LIKE 'Lumba%' AND wavelength > 5728 AND wavelength < 5729""".format(
        element, ion))

# ULB bad for cool stars for 5402.8 line and 5544.6 line ( <4500)
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this star for cool stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (select distinct on (cname) cname,
        feh, teff FROM node_results) n ON (l.cname = n.cname AND 
        TRIM(l.element) = '{0}' AND l.ion = {1} AND 
        l.node LIKE 'ULB%' AND n.teff < 4500 AND 
        (
            (l.wavelength > 5544 AND l.wavelength < 5545) OR
            (l.wavelength > 5402 AND l.wavelength < 5403)
        ))""".format(
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
fig.savefig("figures/{0}{1}/scaled-differential-line-abundances.png".format(
    element.upper(), ion))

fig = ges.plot.solar_comparison(element, ion, scaled=True)
fig.savefig("figures/{0}{1}/scaled-compare-solar.png".format(element.upper(),
    ion))
