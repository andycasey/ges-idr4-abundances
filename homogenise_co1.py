
""" Produce ensemble Co 1 abundances from the Gaia-ESO Survey iDR4 WG11 data """

__author__ = 'Andy Casey <arc@ast.cam.ac.uk>'


import release

database, element, ion = ("arc", "Co", 1)
ges = release.DataRelease(database)

# Lumba routinely gives abundances of like 30+
flag_id = ges.flags.retrieve_or_create("Abundance is non-physical")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1}
    AND (abundance < 1.5 OR abundance > 6)""".format(element, ion))


# Lumba and EPINARBO gives bad values for EMP stars.
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this line for EMP stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        feh FROM node_results) n ON (l.cname = n.cname AND l.element = '{0}'
        AND l.ion = {1} AND n.feh < -2.5 AND (l.node LIKE 'EPINARBO%' OR
            l.node LIKE 'Lumba%'))""".format(element, ion))

# Lumba gives bad results for 5331.4 line
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1} AND
    node LIKE 'Lumba%' AND wavelength > 5331 AND wavelength < 5332""".format(
        element, ion))

# Lumba gives bad results for the 5352 line for cool stars
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this star for cool stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        teff FROM node_results) n ON (l.element = '{0}' AND l.ion = {1} AND
        l.cname = n.cname AND l.node LIKE 'Lumba%' AND n.teff < 4000 AND
        l.wavelength > 5351.5 AND l.wavelength < 5352.5)""".format(element, ion))

flag_id = ges.flags.retrieve_or_create("Not an ideal line for analysis")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1} AND
    wavelength > 5915 AND wavelength < 6006 AND
    (node LIKE 'Lumba%' OR node LIKE 'ULB%')""".format(element, ion))
    
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1} AND
    wavelength > 6490 AND wavelength < 6491""".format(element, ion))


# ULB
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1} AND
    node LIKE 'ULB%' AND wavelength > 6814""".format(element, ion))

flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this star for cool stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        teff FROM node_results) n ON (l.cname = n.cname AND l.element = '{0}'
        AND l.ion = {1} AND l.node LIKE 'ULB%' AND n.teff < 4000
        AND l.wavelength > 6093 AND l.wavelength < 6094)""".format(
        element, ion))
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        teff FROM node_results) n ON (l.cname = n.cname AND l.element = '{0}'
        AND l.ion = {1} AND l.node LIKE 'ULB%' AND n.teff < 4000
        AND l.wavelength > 5530 AND l.wavelength < 5531)""".format(
        element, ion))
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        teff FROM node_results) n ON (l.cname = n.cname AND l.element = '{0}'
        AND l.ion = {1} AND l.node LIKE 'ULB%' AND n.teff < 4000
        AND ((l.wavelength > 5369 AND l.wavelength < 5370) OR 
            (l.wavelength > 5331 AND l.wavelength < 5332)))""".format(
        element, ion))
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        teff FROM node_results) n ON (l.cname = n.cname AND l.element = '{0}'
        AND l.ion = {1} AND l.node LIKE 'ULB%' AND n.teff < 4000
        AND l.wavelength > 5301 AND l.wavelength < 5302)""".format(
        element, ion))

flag_id = ges.flags.retrieve_or_create("Not an ideal line for analysis")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1} AND
    node LIKE 'ULB%' AND wavelength > 5280 AND wavelength < 5281""".format(element, ion))
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1} AND
    node LIKE 'ULB%' AND wavelength > 5176 AND wavelength < 5177""".format(element, ion))


# Calculate biases and apply them.
species_biases = ges.biases.differential(element, ion)
for node in species_biases:
    for wavelength, (bias, sigma, N) in species_biases[node].items():
        rows = ges.biases.apply_offset(element, ion, node, wavelength, -bias)

# Perform the homogenisation.
ges.homogenise.species(element, ion)

