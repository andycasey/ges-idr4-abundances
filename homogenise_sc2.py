
""" Produce ensemble Sc 2 abundances from the Gaia-ESO Survey iDR4 WG11 data """

__author__ = 'Andy Casey <arc@ast.cam.ac.uk>'


import release

database, element, ion = ("arc", "Sc", 2)
ges = release.DataRelease(database)

# Remove lines bluer than 5030 (these come from MyGISFOS only)
flag_id = ges.flags.retrieve_or_create("Not an ideal line for analysis")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1}
    AND wavelength < 5030""".format(element, ion))

# It appears CAUP is not using HFS or something... they show a very large amount
# of scatter relative to all other nodes 
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1} AND
    node LIKE 'CAUP%'""".format(element, ion))

# EPINARBO does bad for metal poor stars up until wavelengths of 5333. And again
# for lines redder than 5668
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this line for EMP stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        feh FROM node_results) n ON (l.cname = n.cname AND l.element = '{0}'
        AND l.ion = {1} AND n.feh < -2.5 AND l.node LIKE 'EPINARBO%'
        AND l.wavelength < 5333)""".format(element, ion))
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        feh FROM node_results) n ON (l.cname = n.cname AND l.element = '{0}'
        AND l.ion = {1} AND n.feh < -2.5 AND l.node LIKE 'EPINARBO%'
        AND l.wavelength > 5668)""".format(element, ion))

# Lumba does bad for cool stars for lines bluer than 5667
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this star for cool stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        teff FROM node_results) n ON (l.cname = n.cname AND l.element = '{0}'
        AND l.ion = {1} AND n.teff < 4000 AND l.node LIKE 'Lumba%' AND
        l.wavelength < 5667)""".format(element, ion))

# Lumba does bad for EMP stars
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this line for EMP stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        feh FROM node_results) n ON (l.cname = n.cname AND l.element = '{0}'
        AND l.ion = {1} AND n.feh < -2.5 AND l.node LIKE 'Lumba%')""".format(
        element, ion))

# ULB remove lines greater than 4.5 abundance for wavelengths < 5040
flag_id = ges.flags.retrieve_or_create("Abundance is non-physical")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1}
    AND abundance > 4.5 AND node LIKE 'ULB%' AND wavelength < 5040"""\
    .format(element, ion))

# Calculate biases and apply them.
species_biases = ges.biases.differential(element, ion)
for node in species_biases:
    for wavelength, (bias, sigma, N) in species_biases[node].items():
        rows = ges.biases.apply_offset(element, ion, node, wavelength, -bias)

# Perform the homogenisation.
ges.homogenise.species(element, ion)
