
""" Produce ensemble Li 1 abundances from the Gaia-ESO Survey iDR4 WG11 data """

__author__ = 'Andy Casey <arc@ast.cam.ac.uk>'


import release

database, element, ion = ("arc", "Li", 1)
ges = release.DataRelease(database)

# Note: OACT results were extracted from the node results file using the script
# in extract_from_node_results.py

# Set all wavelengths to the same value.
ges.execute("""UPDATE line_abundances SET wavelength = 6707.8 WHERE
    element = '{0}' AND ion = {1}""".format(element, ion))
ges.commit()

# The ULB results are well-offset from the Sun, and are offset from all other
# nodes.
flag_id = ges.flags.retrieve_or_create("Spurious results")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE TRIM(element) = '{0}' AND ion = {1}
    AND node LIKE 'ULB%'""".format(element, ion))

# EPINARBO shows large scatter for EMP stars.
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this line for EMP stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        feh FROM node_results) n ON (l.cname = n.cname AND l.element = '{0}'
        AND l.ion = {1} AND n.feh < -2.5 AND l.node LIKE 'EPINARBO%')""".format(
        element, ion))

# OACT provides Li abundances for stars without metallicity. Suspicious.
flag_id = ges.flags.retrieve_or_create(
    "Star does not have a metallicity, therefore no abundance should be available")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) cname,
        feh, teff, logg FROM node_results) n ON (l.cname = n.cname
        AND l.element = '{0}' AND l.ion = {1} AND l.abundance <> 'NaN'
        and n.feh IS NULL)""".format(element, ion))

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
