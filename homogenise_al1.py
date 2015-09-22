
''' Produce ensemble Al 1 abundances from the Gaia-ESO Survey iDR4 WG11 data '''

__author__ = 'Andy Casey <arc@ast.cam.ac.uk>'


import release
import plot


element, ion = ("Al", 1)
ges = release.DataRelease("arc")


# Flag any lines that should be discarded, from specific nodes?
flag_id = ges.flags.retrieve_or_create("Abundance is non-physical")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = '{1}' AND
        (abundance > 8)""".format(element, ion))

# The Al 1 line at 6696.2 from MyGIsFOS is probably a duplicate of 6696.0
flag_id = ges.flags.retrieve_or_create("Line is probably a duplicate print of Al 1 at 6696.0")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = '{1}'
    AND node ILIKE 'MyGIsFOS%' AND wavelength > 6696.1 AND wavelength < 6696.2
    AND flags = 0""".format(element, ion))

# Calculate biases and apply them.
species_biases = ges.biases.differential(element, ion)
for node in species_biases:
    for wavelength, (bias, sigma, N) in species_biases[node].items():
        rows = ges.biases.apply_offset(element, ion, node, wavelength, -bias)

# Produce some figures.
fig = plot.differential_line_abundances(ges._database, element, ion, scaled=True,
    ignore_flags=False)
fig.savefig("figures/{0}{1}/scaled-differential-line-abundances.png".format(
    element.upper(), ion))

fig = plot.compare_solar(ges._database, element, ion, scaled=True, ignore_flags=False)
fig.savefig("figures/{0}{1}/scaled-compare-solar.png".format(element.upper(),
    ion))

# Perform the homogenisation.
ges.homogenise.species(element, ion)

#result = homogenise.species(db, element, ion)

# Produce comparison figures.

# - node-to-node scatter plot
# - benchmark comparison
# - bensby comparison
# - 

raise a 
