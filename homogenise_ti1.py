
""" Produce ensemble Ti 1 abundances from the Gaia-ESO Survey iDR4 WG11 data """

__author__ = 'Andy Casey <arc@ast.cam.ac.uk>'

import release

database, element, ion = ("arc", "Ti", 1)
ges = release.DataRelease(database)

# CAUP shows large scatter below 5700
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1}
    AND node LIKE 'CAUP%' AND wavelength < 5700""".format(
        element, ion))

# CAUP large scatter 5965.9 and 6258.1
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1}
    AND node LIKE 'CAUP%' AND wavelength > 5965.9 AND wavelength < 5966""".format(
        element, ion))

flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1}
    AND node LIKE 'CAUP%' AND wavelength > 6258 AND wavelength < 6259""".format(
        element, ion))

# EPINARBO EMP below 5514 and redder than 5565
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (select distinct on (cname) cname,
        feh FROM node_results) n ON (l.cname = n.cname AND 
        l.element = '{0}' AND l.ion = {1} AND l.node LIKE 'EPINARBO%' AND
        n.feh < -2.5 AND (
            l.wavelength < 5514 OR l.wavelength > 5565
        ))""".format(element, ion))

# Lumba large scatter for 4915.2 and 4926.2 and 4997.1
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1}
    AND node LIKE 'Lumba%' AND (
        (wavelength > 4915 AND wavelength < 4916) OR
        (wavelength > 4926 AND wavelength < 4927) OR
        (wavelength > 4997 AND wavelength < 4998)
    )""".format(element, ion))

# LUMBA EMPs at 5034.6 and 5145.5 and 521s9.7  and 5295.8
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this line for EMP stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (select distinct on (cname) cname,
        feh FROM node_results) n ON (l.cname = n.cname AND 
        l.element = '{0}' AND l.ion = {1} AND l.node LIKE 'ULB%' AND
        n.feh < -2.5 AND (
            (l.wavelength > 5034 AND l.wavelength < 5035) OR
            (l.wavelength > 5145 AND l.wavelength < 5146) OR
            (l.wavelength > 5219 AND l.wavelength < 5220) OR
            (l.wavelength > 5295 AND l.wavelength < 5296)
        ))""".format(element, ion))

# ULB large scatter for 4913.6 and 4915.2
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1}
    AND node LIKE 'ULB%' AND (
        (wavelength > 4913 AND wavelength < 4914) OR
        (wavelength > 4915 AND wavelength < 4916)
    )""".format(element, ion))

# ULB cool stars for 4997.1 and 5009.6 
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this line for EMP stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (select distinct on (cname) cname,
        teff FROM node_results) n ON (l.cname = n.cname AND 
        l.element = '{0}' AND l.ion = {1} AND l.node LIKE 'ULB%' AND
        n.teff < 4250 AND (
            (l.wavelength > 4997 AND l.wavelength < 4998) OR
            (l.wavelength > 5009 AND l.wavelength < 5010)
        ))""".format(element, ion))

# ULB large scatter for 5040.0 and 5043.6 and 5113.4 and 5145.5 and 5147.5 and 5152.2 and 5219.7 and 5338.3 5366.6 and 5384.6 
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1}
    AND node LIKE 'ULB%' AND (
        (wavelength > 5039.5 AND wavelength < 5040.5) OR
        (wavelength > 5043 AND wavelength < 5044) OR
        (wavelength > 5113 AND wavelength < 5114) OR
        (wavelength > 5145 AND wavelength < 5146) OR
        (wavelength > 5147 AND wavelength < 5148) OR
        (wavelength > 5152 AND wavelength < 5153) OR
        (wavelength > 5219 AND wavelength < 5230) OR
        (wavelength > 5338 AND wavelength < 5339) OR
        (wavelength > 5366 AND wavelength < 5367) OR
        (wavelength > 5384 AND wavelength < 5385)
    )""".format(element, ion))

# ULB cool stars for 6126.2
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen in this line for EMP stars")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances l JOIN (select distinct on (cname) cname,
        teff FROM node_results) n ON (l.cname = n.cname AND 
        l.element = '{0}' AND l.ion = {1} AND l.node LIKE 'ULB%' AND
        n.teff < 4250 AND (
            (l.wavelength > 6126 AND l.wavelength < 6127)
        ))""".format(element, ion))

# ULB large scatter for 6395.5 and 6554.2 and 6556.1
flag_id = ges.flags.retrieve_or_create(
    "Large scatter seen for this line from this node")
num_rows = ges.flags.update([flag_id],
    """SELECT id FROM line_abundances WHERE element = '{0}' AND ion = {1}
    AND node LIKE 'ULB%' AND (
        (wavelength > 6395 AND wavelength < 6396) OR
        (wavelength > 6554 AND wavelength < 6555) OR
        (wavelength > 6556 AND wavelength < 6557)
    )""".format(element, ion))


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
