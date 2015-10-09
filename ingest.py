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


def ingest_node_results(filename, kwds):


    connection = pg.connect(**kwds)
    cursor = connection.cursor()
    
    image = fits.open(filename)

    node = filename.split("_")[-1].split(".")[0]

    print("Preparing from {}".format(filename))

    extension = 1 if "MaxPlanck" in filename else 2


    rows = []
    N = len(image[extension].data)
    print(node, N)

    default_row = {"node": node}
    for i, line in enumerate(image[extension].data.base):
        print("Preparing row {0}/{1} from {2}".format(i, N, node))
        row = {}
        row.update(default_row)
        row.update({ k: v for k, v in \
            zip(image[extension].data.dtype.names, line.tolist()) \
            if isinstance(v, (str, unicode)) or np.isfinite(v) })
        
        # Clean up strings because fuck.
        for k in row:
            if isinstance(row[k], (str, unicode)):
                row[k] = row[k].strip()

        rows.append(row)

    print("Inserting from {}".format(filename))
    for i, row in enumerate(rows):
        print("Inserting row {0}/{1} of {2}".format(i, N, node))
        cursor.execute("INSERT INTO node_results ({0}) VALUES ({1});".format(
            ", ".join(row.keys()),
            ", ".join(["%({})s".format(_) for _ in row.keys()])),
        row)

    cursor.close()
    connection.commit()
    connection.close()

    print("Done")



if __name__ == "__main__":

    import psycopg2 as pg
    from astropy.io import fits
    from glob import glob
    import multiprocessing as mp

    kwds = {
        "dbname": "ges-idr4-wg10",
    }
    threads = 8

    connection = pg.connect(**kwds)
    create_tables(connection)

    filenames = glob("/home/arc/Dropbox/WG10/*/Abundances/*.fits")
    
    pool = mp.Pool(threads)
    for filename in filenames:
        pool.apply_async(ingest_node_results, (filename, kwds))

    pool.close()
    pool.join()


    species = [ # As taken from the FITS table headers.
        ("Li", 1),
        ("C", 1),
        ("N", 2),
        ("O", 2),
        ("Ne", 1),
        ("Na", 1),
        ("Mg", 1),
        ("Al", 1),
        ("Si", 1),
        ("S", 1),
        ("Ca", 1),
        ("Sc", 1),
        ("Ti", 1),
        ("V", 1),
        ("Cr", 1),
        ("Mn", 1),
        ("Fe", 1),
        ("Co", 1),
        ("Ni", 1),
        ("Cu", 1),
        ("Zn", 1),
        ("Sr", 1),
        ("Y", 1),
        ("Zr", 1),
        ("Nb", 1),
        ("Mo", 1),
        ("Ru", 1),
        ("Ba", 2),
        ("La", 2),
        ("Ce", 2),
        ("Pr", 2),
        ("Nd", 2),
        ("Sm", 2),
        ("Eu", 2),
        ("Gd", 2),
        ("Dy", 2),
    ]

    import utils

    # Scale the following abundances from these nodes:
    scale_nodes = ["Lumba", "Nice"] # for Mg 1

    cursor = connection.cursor()
    for node in scale_nodes:

        for element, ion in species:
            cursor.execute("""UPDATE node_results set {0}{1} = {0}{1} + {2} WHERE
                TRIM(node) = '{3}'""".format(element.lower(), ion,
                    utils.solar_abundance(element), node))

    cursor.close()
    connection.commit()
    connection.close()
    print("Ingestion complete.")
