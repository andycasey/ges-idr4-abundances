
""" Produce ensemble Nd 2 abundances from the Gaia-ESO Survey iDR4 WG11 data """

__author__ = 'Andy Casey <arc@ast.cam.ac.uk>'

import release

database, element, ion = ("arc", "Nd", 2)
ges = release.DataRelease(database)

# Lumba reports some unusual measurements.
flag_id = ges.flags.retrieve_or_create("Abundance is non-physical")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = '{1}' AND
        (abundance > 10)""".format(element, ion))

# Lumba shows a lot of scatter for the 5842 line (particularly for EMP stars)
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this line for EMP stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = '{1}' AND
        node like 'Lumba%' AND wavelength > 5842 AND wavelength < 5843""".format(element, ion))

# Lumba and ULB disagree a lot for 5276.9 line and 5740.9 line and 4998.5 line
flag_id = ges.flags.retrieve_or_create(
    "Not an ideal line for analysis")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = '{1}' AND
        (wavelength > 4998 AND wavelength < 4999) OR
        (wavelength > 5726 AND wavelength < 5727) OR
        (wavelength > 5740 AND wavelength < 5741)""".format(element, ion))

# Solar biases available for most nodes. Differential biases available for 1/2
# of the remaining lines.

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
