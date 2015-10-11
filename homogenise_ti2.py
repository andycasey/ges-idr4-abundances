
""" Produce ensemble Ti 2 abundances from the Gaia-ESO Survey iDR4 WG11 data """

__author__ = 'Andy Casey <arc@ast.cam.ac.uk>'

import release

database, element, ion = ("arc", "Ti", 2)
ges = release.DataRelease(database)

# CAUP large scatter at 4911.2
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1}
    AND node LIKE 'CAUP%' AND wavelength > 4911 AND wavelength < 4912""".format(
        element, ion))

# CAUP large scatter at 5381.0
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1}
    AND node LIKE 'CAUP%' AND wavelength > 5380.5 AND wavelength < 5381.5""".format(
        element, ion))

# EPINARBO EMPS everything
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this line for cool stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (select distinct on (cname) cname,
        feh FROM node_results) n ON (l.cname = n.cname AND 
        l.element = '{0}' AND l.ion = {1} AND l.node LIKE 'EPINARBO%' AND
        n.feh < -2.5)""".format(element, ion))

# Lumba EMP at 5005.2 at 5381.0 
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this line for EMP stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (select distinct on (cname) cname,
        feh FROM node_results) n ON (l.cname = n.cname AND 
        l.element = '{0}' AND l.ion = {1} AND l.node LIKE 'Lumba%' AND
        n.feh < -2.5 AND (
            (l.wavelength > 5005 AND l.wavelength < 5006) OR
            (l.wavelength > 5380.5 AND l.wavelength < 5381.5)
        ))""".format(element, ion))

# ULB 4798.5 not ideal
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1}
    AND node LIKE 'ULB%' AND wavelength > 4798 AND wavelength < 4799""".format(
        element, ion))

# ULB large scatter 5005.2 and 5129.2 and 5185.9
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1}
    AND node LIKE 'ULB%' AND (
        (wavelength > 5005 AND wavelength < 5006) OR
        (wavelength > 5129 AND wavelength < 5130) OR
        (wavelength > 5185 AND wavelength < 5186)
    )""".format(element, ion))

# ULB cool stars for 5336.8 and 5418.8
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this line for cool stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (select distinct on (cname) cname,
        teff FROM node_results) n ON (l.cname = n.cname AND 
        l.element = '{0}' AND l.ion = {1} AND l.node LIKE 'ULB%' AND
        n.teff < 4000 AND (
            (l.wavelength > 5336 AND l.wavelength < 5337) OR
            (l.wavelength > 5418 AND l.wavelength < 5419)
        ))""".format(element, ion))

# Calculate biases and apply them.
species_biases = ges.biases.differential(element, ion)
for node in species_biases:
    for wavelength, (bias, sigma, N) in species_biases[node].items():
        rows = ges.biases.apply_offset(element, ion, node, wavelength, -bias)

# Perform the homogenisation.
result = ges.homogenise.species(element, ion)


