
""" Produce ensemble Cr 1 abundances from the Gaia-ESO Survey iDR4 WG11 data """

__author__ = 'Andy Casey <arc@ast.cam.ac.uk>'

import release

database, element, ion = ("arc", "Cr", 1)
ges = release.DataRelease(database)

flag_id = ges.flags.retrieve_or_create(
    "Abundance is non-physical")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1}
    AND (abundance < 2 OR abundance > 8)""".format(element, ion))

# CAUP is spurious for 5247.6 and 5348.3
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1}
    AND node LIKE 'CAUP%' AND wavelength > 5247 AND wavelength < 5248""".format(
        element, ion))
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1}
    AND node LIKE 'CAUP%' AND wavelength > 5348 AND wavelength < 5349""".format(
        element, ion))

# Remove EMP lines for all EPINARBO results
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this line for EMP stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (select distinct on (cname) cname,
        feh FROM node_results) n ON (l.cname = n.cname AND 
        l.element = '{0}' AND l.ion = {1} AND l.node LIKE 'EPINARBO%' AND
        n.feh < -2.5)""".format(element, ion))

# Remove ULB for 4789.3 
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1}
    AND node LIKE 'ULB%' AND wavelength > 4789 AND wavelength < 4790""".format(
        element, ion))

# Remove ULB EMP lines for 4936.3 and 5287.2 and 5304.2 and 5318.8 
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this line for EMP stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (select distinct on (cname) cname,
        feh FROM node_results) n ON (l.cname = n.cname AND 
        l.element = '{0}' AND l.ion = {1} AND l.node LIKE 'ULB%' AND
        n.feh < -2.5 AND (
            (l.wavelength > 4936 AND l.wavelength < 4937) OR
            (l.wavelength > 5287 AND l.wavelength < 5288) OR
            (l.wavelength > 5304 AND l.wavelength < 5305) OR
            (l.wavelength > 5318 AND l.wavelength < 5319)
        ))""".format(element, ion))

# Remove ULB fir 5838.7 line and 5844.6
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1}
    AND node LIKE 'ULB%' AND (
        (wavelength > 5838 AND wavelength < 5839) OR
        (wavelength > 5844 AND wavelength < 5845))""".format(element, ion))

# Remove LUMBA for EMP lines below 5204
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this line for EMP stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (select distinct on (cname) cname,
        feh FROM node_results) n ON (l.cname = n.cname AND 
        l.element = '{0}' AND l.ion = {1} AND l.node LIKE 'Lumba%' AND
        n.feh < -2.5 AND l.wavelength < 5204)""".format(element, ion))

# Remove LUMBA for EMP lines between 5214 and 5287.2 
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (select distinct on (cname) cname,
        feh FROM node_results) n ON (l.cname = n.cname AND 
        l.element = '{0}' AND l.ion = {1} AND l.node LIKE 'Lumba%' AND
        n.feh < -2.5 AND l.wavelength > 5214 and l.wavelength < 5287)""".format(element, ion))

# Remove LUMBA EMP lines from 5300 to 5345.0
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (select distinct on (cname) cname,
        feh FROM node_results) n ON (l.cname = n.cname AND 
        l.element = '{0}' AND l.ion = {1} AND l.node LIKE 'Lumba%' AND
        n.feh < -2.5 AND l.wavelength > 5300 and l.wavelength < 5345)""".format(element, ion))

# Remove LUMBA for 5480.5
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1}
    AND node LIKE 'Lumba%' AND (wavelength > 5480 AND wavelength < 5481)""".format(element, ion))


# Remove LUMBA EMP for 5628.6
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this line for EMP stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (select distinct on (cname) cname,
        feh FROM node_results) n ON (l.cname = n.cname AND 
        l.element = '{0}' AND l.ion = {1} AND l.node LIKE 'Lumba%' AND
        n.feh < -2.5 AND l.wavelength > 5628 AND l.wavelength < 5629)""".format(element, ion))

# Remove LUMBA EMP for 5648.3 and 5702.3 
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this line for EMP stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (select distinct on (cname) cname,
        feh FROM node_results) n ON (l.cname = n.cname AND 
        l.element = '{0}' AND l.ion = {1} AND l.node LIKE 'Lumba%' AND
        n.feh < -2.5 AND (
            (l.wavelength > 5648 AND l.wavelength < 5649) OR
            (l.wavelength > 5702 AND l.wavelength < 5703)
        ))""".format(element, ion))

# Remove LUMBA for 5838.7 and 5844.6
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1}
    AND node LIKE 'Lumba%' AND (wavelength > 5838 AND wavelength < 5839)""".format(element, ion))
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1}
    AND node LIKE 'Lumba%' AND (wavelength > 5844 AND wavelength < 5845)""".format(element, ion))



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
