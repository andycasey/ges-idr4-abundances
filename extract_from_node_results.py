
""" Extract OACT results for Li and place them in the line_abundances table """

__author__ = 'Andy Casey <arc@ast.cam.ac.uk>'

import logging
import numpy as np

import release

logger = logging.getLogger("ges")

# We may need to do this many times...
database, remove_existing = ("arc", True)

element, ion, wavelength = ("Li", 1, 6707.8)
node, measurement_type, code = ("OACT", "S", "Unknown")


ges = release.DataRelease(database)

# Remove any existing rows in line_abundances for this element from this node?
if remove_existing:
    logger.info("Deleting existing {0} {1} line abundances from {2}".format(
        element, ion, node))
    ges.execute("""DELETE FROM line_abundances WHERE TRIM(node) = %s
        AND TRIM(element) = %s AND ion = %s""", (node, element, ion))
    ges.commit()

# Get all the details from node results.
node_results = ges.retrieve_table("""SELECT * FROM node_results
    WHERE TRIM(node) = %s""", (node, ))
N = len(node_results)
for i, node_result_row in enumerate(node_results):

    # Create the spectrum_filename_stub
    filenames = node_result_row["filename"].strip().split("|")
    spectrum_filename_stub = ("".join([filenames[0][j] \
        for j in range(max(map(len, filenames))) \
        if len(set([item[j] for item in filenames])) == 1]))[:-5]

    # Create the new row of data.
    #li1 | upper_combined_li1 | e_li1 | nn_li1 | enn_li1 | nl_li1
    upper_column = "upper_{0}{1}".format(element.lower(), ion)
    if upper_column not in node_result_row.dtype.names:
        upper_column = "upper_combined_{0}{1}".format(element.lower(), ion)

    line_abundance_row = {
        "abundance_filename": "GES_iDR4_WG11_{0}.fits".format(node),
        "spectrum_filename_stub": spectrum_filename_stub,
        "node": node,
        "cname": node_result_row["cname"],
        "code": code,
        "object": node_result_row["object"],
        "element": element,
        "ion": ion,
        "wavelength": wavelength,
        "ew": np.nan,
        "e_ew": np.nan,
        "upper_ew": 0,
        "abundance": node_result_row["{0}{1}".format(element.lower(), ion)],
        "e_abundance": node_result_row["e_{0}{1}".format(element.lower(), ion)],
        "upper_abundance": node_result_row[upper_column],
        "measurement_type": measurement_type,
    }
    line_abundance_row["scaled_abundance"] = line_abundance_row["abundance"]

    logger.debug("Inserting row {0}/{1} {2}".format(i + 1, N,
        line_abundance_row.items()))

    ges.execute("""INSERT INTO line_abundances({0}) VALUES ({1})""".format(
        ", ".join(line_abundance_row.keys()),
        ", ".join(["%({})s".format(_) for _ in line_abundance_row.keys()])),
        line_abundance_row)

ges.commit()
logger.info("Done")
