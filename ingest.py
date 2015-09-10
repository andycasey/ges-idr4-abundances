#!/usr/bin/python

""" Ingest all of the node results and put them into a SQLite database. """


import logging
import numpy as np
import os

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("ingest")


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

    rows = []
    metadata = {
        "FILENAME": os.path.basename(filename),
        "NODE": os.path.basename(filename).split("_")[0]
    }
    with open(filename, "r") as fp:
        for i, line in enumerate(fp.readlines()):
            if line.startswith("#"):
                if "CNAME" in line:
                    metadata["CNAME"] = line.split("CNAME")[1].strip()

                elif "CODE" in line:
                    metadata["CODE"] = line.split("CODE")[1].strip()
                    
                elif "OBJECT" in line:
                    metadata["OBJECT"] = line.split("OBJECT")[1].strip()

            else:
                assert len(metadata) == 5

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
                if metadata["NODE"] == "MyGIsFOS":
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
                elif metadata["NODE"] == "EPINARBO":
                    _, element, ion, ew, e_ew, upper_ew, abundance, e_abundance,\
                        upper_abundance, measurement_type = line.split()

                    row.update({
                        "wavelength": wavelength,
                        "element": element,
                        "ion": int(ion),
                        "ew": safe_float(ew),
                        "e_ew": safe_float(ew),
                        "upper_ew": int(upper_ew),
                        "abundance": safe_float(abundance),
                        "e_abundance": safe_float(e_abundance),
                        "upper_abundance": int(upper_abundance),
                        "measurement_type": measurement_type.strip()
                    })

                else:
                    raise WTFError

                rows.append(row)
                logger.debug(row)

    return rows
               

if __name__ == "__main__":

    from glob import glob
    files = glob("data/MyGIsFOS/*.dat") + glob("data/EPINARBO/*.dat")

    for filename in files:
        print(filename)
        foo = parse_line_abundances(filename)
        raise a