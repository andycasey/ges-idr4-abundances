
""" Produce ensemble Fe 2 abundances from the Gaia-ESO Survey iDR4 WG11 data """

__author__ = 'Andy Casey <arc@ast.cam.ac.uk>'

import release

database, element, ion = ("arc", "Fe", 2)
ges = release.DataRelease(database)

# Remove any values below 4 or above 12
flag_id = ges.flags.retrieve_or_create(
    "Abundance is non-physical")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1}
    AND (abundance < 4 OR abundance > 12)""".format(element, ion))

# CAUP shows large scatter at 6516.1
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1}
    AND node LIKE 'CAUP%' AND wavelength > 6156 AND wavelength < 6157""".format(
        element, ion))

# CAUP shows large scatter at wavelengths bluer than 5800. Remove them
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1}
    AND node LIKE 'CAUP%' AND wavelength < 5800""".format(element, ion))


# EPINARBO EMP at
#4893.8
#4993.4
#5000.7
#5132.7
#5197.6
#5264.8
#5284.1
#5325.6
#5337.7
#5414.1
#5425.2
#5525.1
#5534.8
#5991.4
#6084.1
#6113.3
#6149.2
#6238.4
#6247.6
#6416.9
#6432.7
#6456.4
#6516.1

num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (select distinct on (cname) cname,
        feh FROM node_results) n ON (l.cname = n.cname AND 
        l.element = '{0}' AND l.ion = {1} AND l.node LIKE 'EPINARBO%' AND
        n.feh < -2.5 AND (
            (l.wavelength 4893 > AND l.wavelength < 4894) OR
            (l.wavelength 4993 > AND l.wavelength < 4994) OR
            (l.wavelength 5000 > AND l.wavelength < 5001) OR
            (l.wavelength 5132 > AND l.wavelength < 5133) OR
            (l.wavelength 5197 > AND l.wavelength < 5198) OR
            (l.wavelength 5264 > AND l.wavelength < 5265) OR
            (l.wavelength 5284 > AND l.wavelength < 5285) OR
            (l.wavelength 5325 > AND l.wavelength < 5326) OR
            (l.wavelength 5337 > AND l.wavelength < 5338) OR
            (l.wavelength 5414 > AND l.wavelength < 5415) OR
            (l.wavelength 5425 > AND l.wavelength < 5426) OR
            (l.wavelength 5525 > AND l.wavelength < 5526) OR
            (l.wavelength 5534 > AND l.wavelength < 5535) OR
            (l.wavelength 5991 > AND l.wavelength < 5992) OR
            (l.wavelength 6084 > AND l.wavelength < 6085) OR
            (l.wavelength 6113 > AND l.wavelength < 6114) OR
            (l.wavelength 6149 > AND l.wavelength < 6150) OR
            (l.wavelength 6238 > AND l.wavelength < 6239) OR
            (l.wavelength 6247 > AND l.wavelength < 6248) OR
            (l.wavelength 6416 > AND l.wavelength < 6417) OR
            (l.wavelength 6432 > AND l.wavelength < 6433) OR
            (l.wavelength 6456 > AND l.wavelength < 6457) OR
            (l.wavelength 6516 > AND l.wavelength < 6517)
        ))""".format(element, ion))

# Remove MyGISFOS cool stars at 6456.4 (<4250)
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this line for EMP stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (select distinct on (cname) cname,
        teff FROM node_results) n ON (l.cname = n.cname AND 
        l.element = '{0}' AND l.ion = {1} AND l.node LIKE 'My%' AND
        n.teff < 4250 AND (
            l.wavelength > 6456 AND l.wavelength < 6457
        ))""".format(element, ion))

# Remove ULB at 4923.9
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1}
    AND node LIKE 'ULB%' AND wavelength > 4923 AND wavelength < 4924""".format(element, ion))

# Remove ULB at 4993.4
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1}
    AND node LIKE 'ULB%' AND wavelength > 4993 AND wavelength < 4994""".format(element, ion))

# Large scatter for ULB in 5234.6, 5316.6, 5337.7
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1}
    AND node LIKE 'ULB%' AND (
        (wavelength > 5234 AND wavelength < 5235) OR
        (wavelength > 5316 AND wavelength < 5317) OR
        (wavelength > 5337 AND wavelength < 5338)
    )""".format(element, ion))

# Remove ULB cool stars for 5534.8 (<4250)
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this line for cool stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (select distinct on (cname) cname,
        teff FROM node_results) n ON (l.cname = n.cname AND 
        l.element = '{0}' AND l.ion = {1} AND l.node LIKE 'ULB%' AND
        n.teff < 4250 AND (
            l.wavelength > 5334 AND l.wavelength < 5335
        ))""".format(element, ion))

# Remove ULB cool stars for 6149.2 
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this line for cool stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (select distinct on (cname) cname,
        teff FROM node_results) n ON (l.cname = n.cname AND 
        l.element = '{0}' AND l.ion = {1} AND l.node LIKE 'ULB%' AND
        n.teff < 4250 AND (
            l.wavelength > 6149 AND l.wavelength < 6150
        ))""".format(element, ion))

# Remove line at: (not ideal for analysis)
flag_id = ges.flags.retrieve_or_create(
    "Not an ideal line for analysis")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1}
    AND wavelength > 5168.5 AND wavelength < 5169.5""".format(element, ion))


# Calculate biases and apply them.
species_biases = ges.biases.differential(element, ion)
for node in species_biases:
    for wavelength, (bias, sigma, N) in species_biases[node].items():
        rows = ges.biases.apply_offset(element, ion, node, wavelength, -bias)

# Perform the homogenisation.
result = ges.homogenise.species(element, ion)


# Produce some figures.
fig = ges.plot.differential_line_abundances(element, ion, scaled=True)
fig.savefig("figures/{0}{1}/scaled-differential-line-abundances.png".format(
    element.upper(), ion))

fig = ges.plot.solar_comparison(element, ion, scaled=True)
fig.savefig("figures/{0}{1}/scaled-compare-solar.png".format(element.upper(),
    ion))
