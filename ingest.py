#!/usr/bin/python

""" Ingest all of the WG10 node results and put them into a PostgreSQL database. """


import logging
import numpy as np
import os
from astropy.table import Table

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("ges")


def create_tables(connection):

    cursor = connection.cursor()
    with open("schema.sql", "r") as fp:
        cursor.execute(fp.read())
    cursor.close()
    connection.commit()

    return None



if __name__ == "__main__":

    import psycopg2 as pg
    from astropy.io import fits
    from glob import glob

    kwds = {
        "dbname": "ges-idr4-wg10",
    }

    connection = pg.connect(**kwds)
    create_tables(connection)

    filenames = glob("/home/arc/Dropbox/WG10/*/Abundances/*.fits")
    cursor = connection.cursor()
    for filename in filenames:
        image = fits.open(filename)

        print("Preparing from {}".format(filename))

        rows = []
        N = len(image[2].data)
        default_row = {"node": filename.split("_")[-1].split(".")[0]}
        for i, line in enumerate(image[2].data.base):
            print("Preparing row {0}/{1}".format(i, N))
            row = {}
            row.update(default_row)
            row.update({ k: v for k, v in \
                zip(image[2].data.dtype.names, line.tolist()) \
                if isinstance(v, (str, unicode)) or np.isfinite(v) })
            
            # Clean up strings because fuck.
            for k in row:
                if isinstance(row[k], (str, unicode)):
                    row[k] = row[k].strip()

            rows.append(row)

        print("Inserting from {}".format(filename))
        for i, row in enumerate(rows):
            print("Inserting row {0}/{1} of {2}".format(i, N, filename))
            cursor.execute("INSERT INTO node_results ({0}) VALUES ({1});".format(
                ", ".join(row.keys()),
                ", ".join(["%({})s".format(_) for _ in row.keys()])),
            row)

        print("Done")


    # CREATE INDICES
    cursor.close()
    connection.commit()

    print("Ingestion complete.")
