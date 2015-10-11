
""" Produce ensemble Mg 1 abundances from the Gaia-ESO Survey iDR4 WG11 data """

__author__ = 'Andy Casey <arc@ast.cam.ac.uk>'

import release

database, element, ion = ("arc", "Mg", 1)
ges = release.DataRelease(database)

# Remove ULB for 5167.3 line
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1}
    AND node LIKE 'ULB%' AND wavelength > 5167 AND wavelength < 5168""".format(
        element, ion))

# Remove ULB for cool (<4500) or hot (>6000) stars for 5172.7 and 5183.6 line
# and 5528.4 line
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this line for EMP stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (select distinct on (cname) cname,
        teff FROM node_results) n ON (l.cname = n.cname AND 
        l.element = '{0}' AND l.ion = {1} AND l.node LIKE 'ULB%' AND
        (n.teff < 4000 OR n.teff > 6000) AND (
            (l.wavelength > 5172 AND l.wavelength < 5173) OR
            (l.wavelength > 5183 AND l.wavelength < 5184) OR
            (l.wavelength > 5528 AND l.wavelength < 5529)
            ))""".format(
    element, ion))

# Remove ULB for EMP stars for 6318.7 and 6319.2 line
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this line for EMP stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (select distinct on (cname) cname,
        feh FROM node_results) n ON (l.cname = n.cname AND 
        l.element = '{0}' AND l.ion = {1} AND l.node LIKE 'ULB%' AND
        n.feh < -2.5 AND (l.wavelength > 6138 AND l.wavelength < 6139.5))""".format(
    element, ion))

# Remove cool stars for Lumba for all lines
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this line for cool stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (select distinct on (cname) cname,
        teff FROM node_results) n ON (l.cname = n.cname AND 
        l.element = '{0}' AND l.ion = {1} AND l.node LIKE 'Lumba%' AND
        n.teff < 4000)""".format(
    element, ion))

# Remove 5528.4 line for EPINARBO
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1}
    AND node LIKE 'EPINARBO%' AND wavelength > 5528 AND wavelength < 5529""".format(
        element, ion))

# Remobe EMP stars for 6318.7 line for EPINARBO and 6319.2 line 
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this line for EMP stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (select distinct on (cname) cname,
        feh FROM node_results) n ON (l.cname = n.cname AND 
        l.element = '{0}' AND l.ion = {1} AND l.node LIKE 'EPINARBO%' AND
        n.feh < -2.5 AND (l.wavelength > 6138 AND l.wavelength < 6139.5))""".format(
    element, ion))

# Remove 5711.1 line from CAUP
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1}
    AND node LIKE 'CAUP%' AND wavelength > 5711 AND wavelength < 5712""".format(
        element, ion))

# Calculate biases and apply them.
species_biases = ges.biases.differential(element, ion)
for node in species_biases:
    for wavelength, (bias, sigma, N) in species_biases[node].items():
        rows = ges.biases.apply_offset(element, ion, node, wavelength, -bias)

# Perform the homogenisation.
result = ges.homogenise.species(element, ion)

