

import logging
import numpy as np
from astropy.io import fits
import psycopg2 as pg

logging.basicConfig(format='%(asctime)-15s %(level)-8s %(message)s')
logger = logging.getLogger("ges")


def ingest_wg_abundances(db, filename, extension=1, wg=None):

    image = fits.open(filename)
    kwds = {
        "wg": wg or int(filename.split("/")[-1].split("_")[2].strip("WG")),
    }

    # What species are in the table?
    species = [(column[3:-1].title(), int(column[-1])) \
        for column in image[extension].data.dtype.names \
        if column.startswith("NL_") and column.count("_") == 1]
    columns, _,__ = db.execute("SELECT * FROM wg_abundances LIMIT 1")
    
    ingest_rows = []
    N = len(image[extension].data.base)
    for i, row in enumerate(image[extension].data.base):
        print(i, N)

        ingest_row = { k.lower(): v \
            for k, v in zip(image[extension].data.dtype.names, row.tolist()) \
            if k.lower() in columns and (
                isinstance(v, (str, unicode)) or np.isfinite(v)) }
        ingest_row.update(kwds)

        for (element, ion) in species:
            column = "{0}{1}".format(element.upper(), ion)
            upper_column = "UPPER_{0}".format(column)
            if upper_column not in row.dtype.names:
                upper_column = "UPPER_COMBINED_{0}".format(column)
            
            ingest_species_row = ingest_row.copy()
            ingest_species_row.update({
                "element": element,
                "ion": ion,
                "abundance": float(row[column]),
                "e_abundance": float(row["E_" + column]),
                "upper_abundance": int(row[upper_column]),
                "nn": int(row["NN_" + column]),
                "nl": int(row["NL_" + column]),
                "enn": float(row["ENN_" + column])
            })
            if np.isfinite(ingest_species_row["abundance"]):
                ingest_rows.append(ingest_species_row)

    N = len(ingest_rows)
    for i, row in enumerate(ingest_rows):
        print("Ingesting {0}/{1}".format(i, N))
        try:
            db.execute("INSERT INTO wg_abundances ({0}) VALUES ({1})"\
                .format(", ".join(row.keys()), ", ".join(["%({})s".format(_) \
                    for _ in row.keys()])), row)

        except pg.IntegrityError:
            logger.exception("Data integrity error for node-level ingestion:")
            db.execute("rollback")

        else:
            db.commit()


    db.commit()
    image.close()

    return N



def ingest_node_spectrum_abundances(db, filename, extension=2, wg=None):

    image = fits.open(filename)

    # Extract the WG from the filename.
    kwds = {
        "wg": wg or int(filename.split("/")[-1].split("_")[2].strip("WG")),
        "node": image[0].header["NODE1"]
    }

    # What species are in the table?
    species = [(column[3:-1].title(), int(column[-1])) \
        for column in image[extension].data.dtype.names \
        if column.startswith("NL_") and column.count("_") == 1]
    columns, _,__ = db.execute("SELECT * FROM node_spectrum_abundances LIMIT 1")
    
    ingest_rows = []
    N = len(image[extension].data.base)
    for i, row in enumerate(image[extension].data.base):
        print(i, N)

        ingest_row = { k.lower(): v \
            for k, v in zip(image[extension].data.dtype.names, row.tolist()) \
            if k.lower() in columns and (
                isinstance(v, (str, unicode)) or np.isfinite(v)) }
        ingest_row.update(kwds)

        # id, element, ion, flags, abundance, e_abundance, upper_abundance
        for (element, ion) in species:
            column = "{0}{1}".format(element.upper(), ion)
            upper_column = "UPPER_{0}".format(column)
            if upper_column not in row.dtype.names:
                upper_column = "UPPER_COMBINED_{0}".format(column)
            
            ingest_species_row = ingest_row.copy()
            ingest_species_row.update({
                "element": element,
                "ion": ion,
                "abundance": float(row[column]),
                "e_abundance": float(row["E_" + column]),
                "upper_abundance": int(row[upper_column]),
                "enn": float(row["ENN_" + column]),
                "nn": int(row["NN_" + column]),
                "nl": int(row["NL_" + column]) \
                    if np.isfinite(row["NL_" + column]) else 0,
            })
            if np.isfinite(ingest_species_row["abundance"]):
                ingest_rows.append(ingest_species_row)

    N = len(ingest_rows)
    for i, row in enumerate(ingest_rows):
        print("Ingesting {0}/{1}".format(i, N))
        try:
            db.execute("INSERT INTO node_spectrum_abundances ({0}) VALUES ({1})"\
                .format(", ".join(row.keys()), ", ".join(["%({})s".format(_) \
                    for _ in row.keys()])), row)

        except pg.IntegrityError:
            logger.exception("Data integrity error for node-level ingestion:")
            db.execute("rollback")

        else:
            db.commit()

    db.commit()
    image.close()

    return N


def ingest_observations(db, filename, extension=2, wg=None):
    """
    Load in the observational information from each recommended WG file.
    """

    image = fits.open(filename)
    kwds = {
        "wg": wg or int(filename.split("/")[-1].split("_")[2].strip("WG")),
    }

    columns, _,__ = db.execute("SELECT * FROM observations LIMIT 1")
    
    ingest_rows = []
    N = len(image[extension].data.base)
    for i, row in enumerate(image[extension].data.base):
        print(i, N)

        ingest_row = { k.lower(): v \
            for k, v in zip(image[extension].data.dtype.names, row.tolist()) \
            if k.lower() in columns }
        ingest_row.update(kwds)
        ingest_rows.append(ingest_row)


    N = len(ingest_rows)
    for i, row in enumerate(ingest_rows):
        print("Ingesting {0}/{1}".format(i, N))
        try:
            db.execute("INSERT INTO observations ({0}) VALUES ({1})"\
                .format(", ".join(row.keys()), ", ".join(["%({})s".format(_) \
                    for _ in row.keys()])), row)

        except pg.IntegrityError:
            logger.exception("Data integrity error for node-level ingestion:")
            db.execute("rollback")

        else:
            db.commit()

    db.commit()
    image.close()

    return N





from glob import glob
import release

ges = release.DataRelease("ges-idr4-wg15", password="password", host="/tmp")

#ingest_wg_abundances(ges, "data/GES_iDR4_WG13_Recommended.fits")

#for filename in glob("data/GES_iDR4_WG10_*.fits"):
#    if "Recommended" in filename: continue
#    ingest_node_spectrum_abundances(ges, filename)

#for filename in glob("data/GES_iDR4_WG??_Recommended.fits"):
#    ingest_wg_abundances(ges, filename)

for filename in glob("data/GES_iDR4_WG??_Recommended.fits"):
    ingest_observations(ges, filename, 1)

#foo = ingest_node_spectrum_abundances(ges, "data/GES_iDR4_WG11_CAUP.fits")

#ges.execute("CREATE INDEX ON wg_abundances (cname, wg)")
#ges.commit()

