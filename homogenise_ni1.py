
""" Produce ensemble Ni 1 abundances from the Gaia-ESO Survey iDR4 WG11 data """

__author__ = 'Andy Casey <arc@ast.cam.ac.uk>'


import release

database, element, ion = ("arc", "Ni", 1)
ges = release.DataRelease(database)

# Remove lines bluer than 4810
flag_id = ges.flags.retrieve_or_create("Not an ideal line for analysis")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1}
    AND wavelength < 4810""".format(element, ion))


# CAUP results bluer than 5847.0 have very large scatter wrt FEH
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1} AND
    node LIKE 'CAUP%' AND wavelength < 5847.0""".format(element, ion))

# Remove CAUP at 6108.1 line and 6176.8 line
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1} AND
    node LIKE 'CAUP%' AND wavelength > 6108 AND wavelength < 6109""".format(element, ion))

# Remove CAUP line at 6767.8 line
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1} AND
    node LIKE 'CAUP%' AND wavelength > 6767 AND wavelength < 6768""".format(element, ion))

# Remove EPINARBO EMP measurements
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this line for EMP stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        feh FROM node_results) n ON (l.cname = n.cname AND l.element = '{0}'
        AND l.ion = {1} AND n.feh < -2.5 AND l.node LIKE 'EPINARBO%')""".format(element, ion))


# Remove Lumba results bluer than 4913.0
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1} AND
    node LIKE 'Lumba%' AND wavelength < 4913""".format(element, ion))

# Remove 4976,7 line from Lumba
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1} AND
    node LIKE 'Lumba%' AND wavelength < 4976 AND wavelength < 4977""".format(element, ion))

# Remove 5347.7 line from lumba
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1} AND
    node LIKE 'Lumba%' AND wavelength < 5347 AND wavelength < 5348""".format(element, ion))

# Remove 5392.3 line from lumba
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1} AND
    node LIKE 'Lumba%' AND wavelength < 5392 AND wavelength < 5393""".format(element, ion))

# Remove 5468.1 line from Lumba
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1} AND
    node LIKE 'Lumba%' AND wavelength < 5468 AND wavelength < 5469""".format(element, ion))

# Remove 5475.4 line from Lumba
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1} AND
    node LIKE 'Lumba%' AND wavelength < 5475 AND wavelength < 5476""".format(element, ion))

# Remove 5476.9 line from Lumba for stars with FEH > -1
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        feh FROM node_results) n ON (l.cname = n.cname AND n.feh > -1 AND
        l.element = '{0}' AND l.ion = {1} AND l.node LIKE 'Lumba%' AND
        l.wavelength > 5476 AND l.wavelength < 5477)""".format(element, ion))

# Remove 6025.8 line from Lumba
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1} AND
    node LIKE 'Lumba%' AND wavelength < 6025 AND wavelength < 6026""".format(element, ion))

# Remove 6134.0 line from Lumba
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1} AND
    node LIKE 'Lumba%' AND wavelength < 6133.5 AND wavelength < 6134.5""".format(element, ion))

# Remove cool stars for 4980.2 line from lumba andd 5010.9 line and 5035.4 and 5081.1 line and 5084.1 line and 5115.4 line and 6086.3 line and 6086.3 line
# and 6130.1 line
# Remove cool from ALL LUMBA
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this line for cool stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        teff FROM node_results) n ON (l.cname = n.cname AND l.element = '{0}'
        AND l.ion = {1} AND n.teff < 4000 AND l.node LIKE 'Lumba%')""".format(element, ion))

# Remove EMP from all lumba
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this line for EMP stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        feh FROM node_results) n ON (l.cname = n.cname AND l.element = '{0}'
        AND l.ion = {1} AND n.feh < -2.5 AND l.node LIKE 'Lumba%')""".format(element, ion))

# Remove EMP from ULB for all lines except 5476.9 line and 4873.4 line and 4976.3 line and 5035.4 line and 5137.1 line
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this line for EMP stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        feh FROM node_results) n ON (l.cname = n.cname AND l.element = '{0}'
        AND l.ion = {1} AND n.feh < -2.5 AND l.node LIKE 'ULB%' AND (
            (l.wavelength > 6643 AND l.wavelength < 6644) OR
            (l.wavelength > 6482 AND l.wavelength < 6483) OR
            (l.wavelength > 6378 AND l.wavelength < 6379) OR
            (l.wavelength > 6327 AND l.wavelength < 6328) OR
            (l.wavelength > 6230 AND l.wavelength < 6231) OR
            (l.wavelength > 6623.5 AND l.wavelength < 6624.5) OR
            (l.wavelength > 6204 AND l.wavelength < 6205) OR
            (l.wavelength > 6108 AND l.wavelength < 6109) OR
            (l.wavelength > 5486.5 AND l.wavelength < 5487.5) OR
            (l.wavelength > 5435 AND l.wavelength < 5436) OR
            (l.wavelength > 5157.5 AND l.wavelength < 5158.5) OR
            (l.wavelength > 4976 AND l.wavelength < 4977) OR
            (l.wavelength > 4811.5 AND l.wavelength < 4812.5)
        ))""".format(element, ion))

# Remove 4976.7 line from ULB
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1} AND
    node LIKE 'ULB%' AND wavelength < 4976 AND wavelength < 4977""".format(element, ion))


# Calculate biases and apply them.
species_biases = ges.biases.differential(element, ion)
for node in species_biases:
    for wavelength, (bias, sigma, N) in species_biases[node].items():
        rows = ges.biases.apply_offset(element, ion, node, wavelength, -bias)

# Perform the homogenisation.
ges.homogenise.species(element, ion)

