
""" Produce ensemble Zn 1 abundances from the Gaia-ESO Survey iDR4 WG11 data """

__author__ = 'Andy Casey <arc@ast.cam.ac.uk>'


import release

database, element, ion = ("arc", "Zn", 1)
ges = release.DataRelease(database)

# ULB returned a large abundance value for one star
flag_id = ges.flags.retrieve_or_create("Abundance is non-physical")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = '{1}' AND
        (abundance > 7.5)""".format(element, ion))

# Only two lines present (4810.5 and 6362.3) and even though the abundances are
# generally offset from each other, they look sensible.

# Calculate biases and apply them.
species_biases = ges.biases.differential(element, ion)
for node in species_biases:
    for wavelength, (bias, sigma, N) in species_biases[node].items():
        rows = ges.biases.apply_offset(element, ion, node, wavelength, -bias)

# Perform the homogenisation.
ges.homogenise.species(element, ion)

# Produce some figures.
fig = plot.differential_line_abundances(ges._database, element, ion, scaled=True,
    ignore_flags=False)
fig.savefig("figures/{0}{1}/scaled-differential-line-abundances.png".format(
    element.upper(), ion))

fig = plot.compare_solar(ges._database, element, ion, scaled=True, ignore_flags=False)
fig.savefig("figures/{0}{1}/scaled-compare-solar.png".format(element.upper(),
    ion))
