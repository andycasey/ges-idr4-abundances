
""" Produce ensemble S 1 abundances from the Gaia-ESO Survey iDR4 WG11 data """

__author__ = 'Andy Casey <arc@ast.cam.ac.uk>'


import release

database, element, ion = ("arc", "S", 1)
ges = release.DataRelease(database)


flag_id = ges.flags.retrieve_or_create("Abundance is non-physical")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE TRIM(element) = '{0}' AND ion = '{1}' AND
        (abundance > 9)""".format(element, ion))

# ULB and EPINARBO show very different results (~0.4 dex offset) for lines bluer
# than 6757.2. Thus, we will just use the 6757.2 line.
flag_id = ges.flags.retrieve_or_create(
    "Not an ideal line for analysis")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE TRIM(element) = '{0}' AND ion = '{1}'
    AND wavelength < 6757 AND abundance != 'NaN'""".format(element, ion))

# ULB for EMP stars
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this line for EMP stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (select distinct on (cname) cname,
        feh FROM node_results) n ON (l.cname = n.cname AND 
        TRIM(l.element) = '{0}' AND l.ion = {1} AND l.node LIKE 'ULB%' AND
        n.feh < -2.5 AND (
            (l.wavelength > 6757 AND l.wavelength < 6758)
        ))""".format(element, ion))

# Calculate biases and apply them.
species_biases = ges.biases.differential(element, ion)
for node in species_biases:
    for wavelength, (bias, sigma, N) in species_biases[node].items():
        rows = ges.biases.apply_offset(element, ion, node, wavelength, -bias)

# Perform the homogenisation.
ges.homogenise.species(element, ion)

