#!/usr/bin/python

""" Ingest all of the node results and put them into a SQLite database. """


import logging
import numpy as np
import os
from astropy.table import Table

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("ingest")

def load_benchmark_lines(filename):

    return Table.read(filename, format="ascii", 
        names=("wavelength", "element",  "ion", "loggf", "chi", "ges_flags",
            "bm_fg_dwarf", "bm_k_dwarf", "bm_mp_flag", "bm_fgk_giant",
            "bm_m_giant"))   




def parse_line_abundances(filename):
    """
    Parse line abundance information from a single filename.

    :param filename:
        The filename that contains individual line abundances from a node.

    :type filename:
        str

    :returns:
        A list of dictionaries that account for each line abundance given.
    """

    def safe_float(s):
        try:
            s = float(s)
        except (TypeError, ValueError):
            return np.nan
        else:
            return s

    def safe_int(s):
        try:
            s = int(s)
        except (TypeError, ValueError):
            return 0
        else:
            return s

    rows = []
    basename = os.path.basename(filename)
    #EPINARBO_uv_11504236+0145494_580.0_16.wg11_abun_idr4.final.dat
    stub = "_".join(basename.split(".wg11")[0].split("_")[1:])
    metadata = {
        "abundance_filename": basename,
        "spectrum_filename_stub": stub,
        "node": os.path.basename(filename).split("_")[0]
    }
    with open(filename, "r") as fp:
        for i, line in enumerate(fp.readlines()):
            if line.startswith("#"):
                if "CNAME" in line:
                    metadata["cname"] = line.split("CNAME")[1].strip()

                elif "CODE" in line:
                    _ = ["CODE", "CODES"]["CODES" in line]    
                    metadata["code"] = line.split(_)[1].strip()
                
                elif "Code" in line:
                    # ULB
                    metadata["code"] = line[6:]


                elif "OBJECT" in line:
                    metadata["object"] = line.split("OBJECT")[1].strip()

            elif line.startswith("LAMBDA"):
                continue

            else:

                assert len(metadata) == 6

                # Split up the row.
                """
                #CNAME 12423899-1305127
                #OBJECT U_1_b_18_33
                #CODE AUTO_EW v1.0
                #
                # (F9.4,2X,A2,2X,I1,2X,F7.2,2X,F7.2,2X,I1,F5.2,2X,F5.2,2X,I1,2X,A2)
                LAMBDA  ELEMENT  ION  EW  e_EW  Upper_EW  ABUND  e_ABUND  Upper_ABUND  MEAS_TYPE
                3131.0667  Be  2  1100.00  1100.00  0  89.00  50.00  0  SS
                """
                
                wavelength = float(line[:9])
                # Convert wavelengths to Angstroms if needed.
                if 1000 > wavelength: wavelength *= 10.

                row = {}
                row.update(metadata)
                if metadata["node"] == "MyGIsFOS":
                    row.update({
                        "wavelength": wavelength,
                        "element": line[11:14].strip(),
                        "ion": int(line[15:16]),
                        "ew": safe_float(line[18:25]),
                        "e_ew": safe_float(line[27:34]),
                        "upper_ew": int(line[36]),
                        "abundance": safe_float(line[37:42]),
                        "e_abundance": safe_float(line[44:49]),
                        "upper_abundance": int(line[51]),
                        "measurement_type": line[53:55].strip()
                    })
                elif metadata["node"] == "EPINARBO":
                    _, element, ion, ew, e_ew, upper_ew, abundance, e_abundance,\
                        upper_abundance, measurement_type = line.split()

                    row.update({
                        "wavelength": wavelength,
                        "element": element,
                        "ion": int(ion),
                        "ew": safe_float(ew),
                        "e_ew": safe_float(e_ew),
                        "upper_ew": int(upper_ew),
                        "abundance": safe_float(abundance),
                        "e_abundance": safe_float(e_abundance),
                        "upper_abundance": int(upper_abundance),
                        "measurement_type": measurement_type.strip()
                    })

                elif metadata["node"] == "Lumba":
                    _, element, ion, ew, e_ew, upper_ew, abundance, e_abundance,\
                        upper_abundance, measurement_type = line.split()

                    row.update({
                        "wavelength": wavelength,
                        "element": element,
                        "ion": int(ion),
                        "ew": safe_float(ew),
                        "e_ew": safe_float(e_ew),
                        "upper_ew": safe_int(upper_ew),
                        "abundance": safe_float(abundance),
                        "e_abundance": safe_float(e_abundance),
                        "upper_abundance": safe_int(upper_abundance),
                        "measurement_type": measurement_type.strip()
                    })

                elif metadata["node"] in ("ULB", "CAUP"):
                    _, element, ion, ew, e_ew, upper_ew, abundance, e_abundance,\
                        upper_abundance, measurement_type = line.split()

                    row.update({
                        "wavelength": wavelength,
                        "element": element,
                        "ion": int(ion[0]),
                        "ew": safe_float(ew),
                        "e_ew": safe_float(e_ew),
                        "upper_ew": safe_int(upper_ew),
                        "abundance": safe_float(abundance),
                        "e_abundance": safe_float(e_abundance),
                        "upper_abundance": safe_int(upper_abundance),
                        "measurement_type": measurement_type.strip()
                    })

                elif metadata["node"] == "Vilnius":
                    _, element, ion, ew, e_ew, upper_ew, abundance, e_abundance,\
                        upper_abundance, measurement_type = line.split()

                    row.update({
                        "wavelength": wavelength,
                        "element": element,
                        "ion": int(ion),
                        "ew": safe_float(ew),
                        "e_ew": safe_float(e_ew),
                        "upper_ew": safe_int(upper_ew),
                        "abundance": safe_float(abundance),
                        "e_abundance": safe_float(e_abundance),
                        "upper_abundance": safe_int(upper_abundance),
                        "measurement_type": measurement_type.strip()
                    })

                else:
                    raise WTFError
                    
                assert 3 > row["ion"]

                # Assume no scaling/homogenisation for any abundances.
                row["scaled_abundance"] = row["abundance"]
                rows.append(row)
                logger.debug(row)

    return rows


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
    #    "host": "/tmp/",
        "dbname": "ges-idr4-wg11-differential-bias",
    #    "user": "arc",
    #    "password": os.environ.get('PSQL_PW', None)
    }

    connection = pg.connect(**kwds)


    create_tables(connection)
    
    files = glob("data/*/*.dat")

    cursor = connection.cursor()
    for filename in files:
        print(filename)
        line_abundances = parse_line_abundances(filename)
        if len(line_abundances) == 0: continue

        k = line_abundances[0].keys()
        cursor.executemany("""INSERT INTO line_abundances({0}) VALUES ({1})"""\
            .format(", ".join(k), ", ".join(["%({})s".format(_) for _ in k])),
            line_abundances)

    cursor.execute("create index cname_species_index on line_abundances (cname, element, ion);")
    cursor.execute("create index wavelength_index on line_abundances (wavelength);")
    cursor.close()
    connection.commit()
    
    filenames = glob("data/*.fits")
    cursor = connection.cursor()
    for filename in filenames:
        image = fits.open(filename)

        print("Preparing from {}".format(filename))

        rows = []
        default_row = {"node": filename.split("_")[-1].split(".")[0]}
        for i, line in enumerate(image[2].data.base):
            row = {}
            row.update(default_row)
            row.update({ k: v for k, v in zip(image[2].data.dtype.names, line.tolist()) if isinstance(v, (str, unicode)) or np.isfinite(v) })
            # Clean up strings because fuck.
            for k in row:
                if isinstance(row[k], (str, unicode)):
                    row[k] = row[k].strip()

                # EPINARBO:
                if k.startswith("ENN_") and (isinstance(row[k], (str, unicode))\
                    or not np.isfinite(row[k]) or row[k] < 0):
                    row[k] = 0

            if "STAR" in row:
                del row["STAR"]
            rows.append(row)

        print("Inserting from {}".format(filename))

        for row in rows:
            cursor.execute("INSERT INTO node_results ({0}) VALUES ({1});".format(
                ", ".join(row.keys()), ", ".join(["%({})s".format(_) for _ in row.keys()])),
            row)

        print("Done")


    # CREATE INDICES

    cursor.close()
    connection.commit()

    '''
    # Input data from the GES line list.
    benchmark_lines = load_benchmark_lines("benchmark_lines.dat")

    cursor = connection.cursor()
    for line in benchmark_lines:
        print("Inserting info for line {}".format(line))
        # Query for indices matching that.
        data = cursor.execute("""SELECT id FROM line_abundances WHERE AND element = %s
            AND ion = %s AND wavelength > %s AND wavelength < %s""",
            (line["element"], line["ion"], line["wavelength"] - 0.1,
                line["wavelength"] + 0.1))
        ids = cursor.fetchall()
        if len(ids) > 0:
            ids = tuple([int(_[0]) for _ in ids])
            cursor.execute("""UPDATE line_abundances SET loggf = %s, chi = %s
                WHERE id in %s""", (line["loggf"], line["chi"], ids))

    cursor.close()
    connection.commit()
    '''

    print("Ingestion complete.")
