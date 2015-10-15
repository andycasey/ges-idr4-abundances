

""" Script to remove spurious WG12 recommended measurements. """

import logging
import psycopg2 as pg
import numpy as np
from collections import Counter
from astropy.io import fits

import release


logging.basicConfig(format='%(asctime)-15s %(level)-8s %(message)s')
logger = logging.getLogger("ges")


ges = release.DataRelease("ges-idr4-wg15", password="password", host="/tmp")


# Load in the template of CNAMEs.
wg15_template = fits.open("GES_iDR4_WG15_Recommended_FFFT_30072015.fits")

# Exclude any duplicate CNAMEs
mask = np.ones(len(wg15_template[1].data), dtype=bool)
for duplicate_cname in (k for k, v in \
    Counter(wg15_template[1].data["CNAME"]).items() if v > 1):
    mask[np.where(wg15_template[1].data["CNAME"] == duplicate_cname)[0][1:]] \
        = False

wg15_template[1].data = wg15_template[1].data[mask]


# Flag abundances from WG12 for particular elements, etc. We do this here so
# that these measurements won't be considered at all, so that abundances from
# other WGs can take precedence over these measurements.
N = ges.update("""UPDATE wg_abundances SET flags = 1
    WHERE wg = 12 AND (
        (element = 'Mn' AND ion = '1') OR
        (element = 'Co' AND ion = '1') OR
        (element = 'Si' AND ion = '2') OR
        (trim(element) = 'V' AND ion = '1')
    )""")
assert N > 0

N = ges.update("""UPDATE wg_abundances SET flags = 1
    WHERE wg = 12 AND element = 'La' AND ion = '2'
        AND cname in (
            SELECT DISTINCT ON (cname) cname
            FROM observations WHERE logg > 3)""")
assert N > 0

# We also need to do WG10 setup-to-setup corrections here because WG10 is never
# competent enough to do their own homogenisation properly.

# WG10 - WG11 biases and corrections.

# These are HR10|HR21 corrections to put on the WG11 scale.
# Use a sanity check to ensure we do not supply these corrections multiple times.

# Uncorrected:
# filename = gar_09595100+2955171_H548.8.fit|gar_09595100+2955171_H875.7.fit
# abundance of Al 1 in HR10|HR21: 6.36752462387085
expected_abundance = 6.36752462387085
wg10_abundance = ges.retrieve("""SELECT abundance FROM wg_abundances
    WHERE element = 'Al' AND 
        wg = 10 AND
        filename = 'gar_09595100+2955171_H548.8.fit|gar_09595100+2955171_H875.7.fit'
        """)[0][0]

assert float(wg10_abundance) == expected_abundance


wg10_hr1021_to_wg11 = {
    # m * [Fe/H] + b
    ("Mg", 1): (0, 0.0625029),
    ("Al", 1): (0, -0.0765719),
    ("Si", 1): (0, -0.0911155),
    ("Ca", 1): (0, -0.00802088),
    ("Ca", 2): (0.180948, -0.169038),
    ("Ti", 1): (0, 0.0101333),
    ("Ti", 2): (0, -0.0392075),
    ("Cr", 1): (-0.283576, -0.00421249),
    ("Mn", 1): (-0.463962, -0.208325),
    ("Co", 1): (0, -0.0690413),
    ("Ni", 1): (0, 0.0482459)
}

for (element, ion), (m, b) in wg10_hr1021_to_wg11.items():
    N = ges.update("""UPDATE wg_abundances
        SET abundance = wg_abundances.abundance + (O.feh * %s + %s)
        FROM wg_abundances W JOIN observations O ON (
            W.filename = O.filename
            )
        WHERE W.wg = 10 AND W.element = %s AND
            W.ion = %s AND O.setup = 'HR10|HR21'""",
            (m, b, element, ion))
    assert N > 0


hr15n_to_wg11 = {
    ("C", 1): -0.262,
    ("Al", 1): -0.37,
    ("Si", 1): -0.163,
    ("Ca", 1): -0.035,
    ("Ti", 1): -0.180,
    ("Co", 1): -0.052,
    ("Ni", 1): -0.035,
    ("Ba", 2): -0.444
}

for (element, ion), offset in hr15n_to_wg11.items():
    N = ges.update("""UPDATE wg_abundances
        SET abundance = abundance + %s
        FROM wg_abundances W JOIN observations O ON (
            W.wg = 10 AND TRIM(W.element) = %s AND W.ion = %s AND
            TRIM(W.filename) = TRIM(O.filename) AND TRIM(O.setup) = 'HR15N')""",
        (offset, element, ion))
    assert N > 0

hr9b_to_wg11 = {
    ("Cr", 1): -0.12,
    ("Ti", 1): 0.02
}
for (element, ion), offset in hr9b_to_wg11.items():
    N = ges.update("""UPDATE wg_abundances
        SET abundance = abundance + %s
        FROM wg_abundances W JOIN observations O ON (
            W.wg = 10 AND TRIM(W.element) = %s AND W.ion = %s AND
            TRIM(W.filename) = TRIM(O.filename) AND TRIM(O.setup) = 'HR9B')""",
        (offset, element, ion))
    assert N > 0


ges.commit()


# Generate a full list of species available (from any star).
species = ges.retrieve_table("""SELECT DISTINCT ON (element, ion) element, ion
    FROM wg_abundances WHERE abundance <> 'NaN' AND flags = 0""")




def insert_survey_abundance(cname, wg, element, ion):

    element = element.strip() # To match properly against the database.

    # These are the columns that will be transferred from the wg-level
    # abundances to the survey-level abundances.
    columns = ("id", "wg", "cname", "element", "ion", "abundance",
        "e_abundance", "enn", "nn", "nl", "upper_abundance", "flags")

    # Here we have the issue of multiple spectra per CNAME for a given WG.
    N_spectra = ges.retrieve("""SELECT COUNT(*)
        FROM (
            SELECT DISTINCT ON (filename) filename FROM wg_abundances
            WHERE cname = %s AND wg = %s AND
                TRIM(element) = %s AND ion = %s
                AND abundance <> 'NaN' AND flags = 0
        ) AS temp""", (cname, wg, element, ion))[0][0]


    # If there is only one spectrum for this CNAME, the problem is trivial.
    if N_spectra == 0:
        # Nothing to do.
        return
    elif N_spectra == 1:

        ges.execute("""INSERT INTO survey_abundances
                (id_provenance, wg_provenance, {common_columns})
            SELECT {columns} FROM wg_abundances AS T """.format(
                element=element, ion=ion, columns=", ".join(columns),
                common_columns=", ".join(columns[2:])) + \
            """WHERE T.cname = %s AND T.wg = %s
                AND TRIM(T.element) = %s AND T.ion = %s AND T.abundance <> 'NaN'
                AND T.flags = 0
            """, (cname, wg, element, ion))
    else:

        # When there are multiple, we should take a weighted average of the
        # values.
        data = ges.retrieve_table("""SELECT * FROM wg_abundances
            WHERE cname = %s AND wg = %s AND TRIM(element) = %s AND ion = %s
            AND abundance <> 'NaN' AND flags = 0""", (cname, wg, element, ion))

        x = data["abundance"].astype(float)
        x_err = data["e_abundance"].astype(float)

        # If all x and x_err are the same, we should actually treat this as if
        # it were N_spectra == 1.
        if np.std(x) == 0 and np.std(x_err) == 0:
            ges.execute("""INSERT INTO survey_abundances
                (id_provenance, wg_provenance, {common_columns})
            SELECT {columns} FROM wg_abundances AS T """.format(
                element=element, ion=ion, columns=", ".join(columns),
                common_columns=", ".join(columns[2:])) + \
            """WHERE T.cname = %s AND T.wg = %s
                AND TRIM(T.element) = %s AND T.ion = %s AND T.abundance <> 'NaN'
                AND T.flags = 0 AND T.filename = %s
            """, (cname, wg, element, ion, data["filename"][0]))

            return True

        # Clean up errors, because we have to.
        x_err[0 >= x_err] = np.nan
        no_errors_reported = False
        if not np.all(np.isfinite(x_err)):
            logger.warn("Non-finite or zero errors reported by WG{0} for "\
                "{1} {2} abundance for {3}".format(wg, element, ion, cname))

            # If they are all non-finite, let's just take mean +/- sigma/sqrt(N)
            if not np.any(np.isfinite(x_err)):
                no_errors_reported = True

            else:
                # If just some are non-finite, fill the rest with average error
                x_err[~np.isfinite(x_err)] = np.mean(x_err[np.isfinite(x_err)])

        assert np.all(np.isfinite(x))
        assert np.all(np.isfinite(x_err)) or no_errors_reported
        assert np.all(x_err > 0) or no_errors_reported
        assert np.all(1 > data["upper_abundance"]) # Future-Andy problem.
    
        if no_errors_reported:
            abundance = np.mean(x)
            error = np.std(x)/np.sqrt(len(x))

        else:
            # When there are multiple measurements, we should take a weighted
            # average of the available data.

            weights = 1.0/(x_err**2)
            weights /= sum(weights)

            abundance = np.sum(weights * x)
            error = np.sqrt(np.sum(weights**2 * x_err**2))

        assert np.isfinite([abundance, error]).all()
        
        ges.execute("""INSERT INTO survey_abundances
            (id_provenance, wg_provenance,
                cname, element, ion, abundance, e_abundance, upper_abundance,
                enn, nn, nl)
            VALUES (-1, %(wg)s, %(cname)s, %(element)s, %(ion)s, %(abundance)s,
                %(e_abundance)s, %(upper_abundance)s, %(enn)s, %(nn)s, %(nl)s)
            """, {
                "wg": wg,
                "cname": cname,
                "element": element,
                "ion": ion,
                "abundance": abundance,
                "e_abundance": error,
                "upper_abundance": 0,
                "enn": data["enn"][0], # TODO: CHECK.
                "nn": np.nanmax(data["nn"]),
                "nl": np.nanmax(data["nl"])
            })

    return True


# For each CNAME, see which WG was used for the stellar parameter determination.
for i, (cname, recommended_wg) in enumerate(zip(
    wg15_template[1].data["CNAME"], wg15_template[1].data["REC_WG"])):
    
    try:
        recommended_wg = int(recommended_wg.strip("WG "))
    except ValueError:
        print("Error in parsing recommended WG: {0}".format(recommended_wg))
    
    # Check if there are *any* abundances.
    abundances = ges.retrieve_table("""SELECT DISTINCT ON (wg, element, ion)
        wg, element, ion FROM wg_abundances
        WHERE cname = %s AND abundance <> 'NaN' AND flags = 0""", (cname, ))

    all_abundances = ges.retrieve_table("""SELECT * FROM wg_abundances
        WHERE cname = %s AND abundance <> 'NaN' AND flags = 0""", (cname, ))

    # In the simplest case, there are no abundances.
    if abundances is None: continue

    # Were there abundances from the WGs that we did not expect?
    contributing_wgs = sorted(set(abundances["wg"]))
    
    # In the second simplest case, there are abundances only from one WG.
    if len(contributing_wgs) == 1:

        if contributing_wgs[0] != recommended_wg:
            print("Warning: only one WG contributing and it's not as we "\
                "expected: {0} != {1}".format(contributing_wgs[0], recommended_wg))

        # The implicit assumption is that there is only one spectrum for each
        # abundance provided.
        for row in abundances:
            insert_survey_abundance(cname, contributing_wgs[0], row["element"],
                row["ion"])

        # Continue on to the next CNAME, because there is nothing else to do:
        continue

    else:

        # Prioritise by WG, then group measurements together.
        priority = (abundances["wg"] == recommended_wg)
        for row in abundances[priority]:
            insert_survey_abundance(cname, recommended_wg, row["element"],
                row["ion"])


        # Any missing species?
        cname_species = [(r["element"], r["ion"]) for r in abundances]
        missing = set(cname_species).difference([(r["element"], r["ion"]) \
            for r in abundances[priority]])

        missing = missing.difference([("Fe", 1)])
        print("MISSING WG11 Fe 1")

        if missing:

            # For each missing species, is it just one WG? If so, awesome.
            for element, ion in missing:
                    
                mask = (abundances["element"] == element) \
                    * (abundances["ion"] == ion)

                if len(abundances["wg"][mask]) == 1:
                    insert_survey_abundance(cname, abundances["wg"][mask][0],
                        element, ion)

                else:
                    # At this stage we should be going through THE PRIORITIES,
                    # however, it seems in these 43 cases it is only WG10 and
                    # WG11 available. In these cases WG11 always takes priority.

                    assert len(abundances["wg"][mask]) == 2 and \
                        10 in list(abundances["wg"][mask]) and \
                        11 in list(abundances["wg"][mask])

                    insert_survey_abundance(cname, 11, element, ion)


# Apply biases or corrections from different WGs to the *survey_abundances*
# table, based on provenance, etc.

# These corrections are current as of October 15th, and came from:
# https://docs.google.com/document/d/1LHO-_4s3ufgCFIFcQtW0M5xO_LQoiyZwOKjeNKq49x0/edit

# WG13 biases (WG11 - WG13) based on stars in common clusters.
wg13_biases = {
    ("C", 1): [
        (0.05, 2), # NGC 2516
        (0.09, 5), # NGC 2547
        (0.06, 6)  # NGC 6633
        ],
    ("O", 1): [
        (0.08, 2), #    .
        (0.48, 5), #    .
        (0.10, 6)  #    .
    ],
    ("Mg", 1): [
        (0.20, 2),
        (0.35, 5),
        (0.05, 6)
    ],
    ("Sc", 2): [
        (0.09, 2),
        (0.07, 5),
        (0.01, 6)
    ],
    ("Fe", 2): [
        (0.08, 2),
        (0.07, 5),
        (0.13, 6)
    ]
}

for (element, ion), offsets in wg13_biases.items():

    # Calculate average offset, as weighted by number of members.
    offsets = np.array(offsets)
    weights = offsets[:, 1]/offsets[:, 1].sum()
    average_bias = np.sum(offsets[:, 0] * weights)

    N = ges.update("""UPDATE survey_abundances SET abundance = abundance + %s
        WHERE wg_provenance = 13 AND element = %s AND ion = %s""",
        (average_bias, element, ion))
    assert N > 0

    print("Average bias for WG11-WG13 applied for {0} {1}: {2}".format(
        element, ion, average_bias))


# WG12 biases (WG11 - WG12) based on 1339 stars in common between the two.
# We need to compute these offsets.



# WG10 biases were applied on a per-setup basis before survey-level homogenisation
# because they are incompetent.




ges.commit() 


# Set minimum uncertainty on abundances of 0.1 dex for non-WG11 data if the
# S/N is less than 20.
poor_snr_cnames \
    = list(wg15_template[1].data["CNAME"][wg15_template[1].data["SNR"] < 20])
ges.execute("""UPDATE survey_abundances
    SET e_abundance = GREATEST(e_abundance, 0.1)
    WHERE wg_provenance <> 11 AND cname IN %s""", (tuple(poor_snr_cnames), ))

# For all stars, set the uncertainty floor for each element/species to be the
# minimum achievable from WG11.
species_from_other_wgs = ges.retrieve_table("""SELECT DISTINCT ON
        (element, ion) element, ion
    FROM survey_abundances WHERE wg_provenance <> 11""")

for row in species_from_other_wgs:
    element, ion = row["element"], row["ion"]

    # Minimum abundance from WG11?
    min_abundance = ges.retrieve("""SELECT e_abundance FROM survey_abundances
        WHERE wg_provenance = 11 AND element = %s AND ion = %s
        ORDER BY e_abundance ASC LIMIT 1""", (element, ion))

    if len(min_abundance) == 0:
        logger.warn("No minimum abundance for {0} {1}".format(element, ion))
        continue

    min_abundance = float(min_abundance[0][0])
    ges.execute("""UPDATE survey_abundances
        SET e_abundance = GREATEST(e_abundance, %s)
        WHERE wg_provenance <> 11 AND element = %s AND ion = %s""",
        (min_abundance, element, ion))

# Export the data into the template fits file.


# Write the table.


