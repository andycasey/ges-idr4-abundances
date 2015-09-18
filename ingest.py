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
                    
                elif "OBJECT" in line:
                    metadata["object"] = line.split("OBJECT")[1].strip()

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
                        "e_ew": safe_float(ew),
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
                        "e_ew": safe_float(ew),
                        "upper_ew": safe_int(upper_ew),
                        "abundance": safe_float(abundance),
                        "e_abundance": safe_float(e_abundance),
                        "upper_abundance": safe_int(upper_abundance),
                        "measurement_type": measurement_type.strip()
                    })

                else:
                    raise WTFError

                # Assume no scaling/homogenisation for any abundances.
                row["scaled_abundance"] = row["abundance"]
                rows.append(row)
                logger.debug(row)

    return rows


def create_tables(connection):

    cursor = connection.cursor()

    cursor.execute("""DROP TABLE IF EXISTS line_abundances;""")
    cursor.execute("""DROP TABLE IF EXISTS node_results;""")
    cursor.execute("DROP TABLE IF EXISTS line_abundance_flags;")

    cursor.execute("""CREATE TABLE line_abundance_flags(
        description char(120) not null);""")
    cursor.execute("""ALTER TABLE line_abundance_flags ADD COLUMN id BIGSERIAL PRIMARY KEY;""")

    cursor.execute("""CREATE TABLE line_abundances(
        node char(10) not null,
        wavelength real not null,
        element char(2) not null,
        ion integer not null,
        ew real,
        e_ew real,
        upper_ew int,
        abundance real,
        e_abundance real,
        upper_abundance int,
        measurement_type char(2) not null,
        cname char(16) not null,
        code char(30) not null,
        object char(21) not null,
        abundance_filename char(140) not null,
        spectrum_filename_stub char(140) not null,
        scaled_abundance real default 'NaN',
        flags int default 0
        );""")
    cursor.execute("ALTER TABLE line_abundances ADD COLUMN id BIGSERIAL PRIMARY KEY;")

    cursor.execute("""CREATE TABLE node_results(
        node char(10) not null,
        cname char(16) not null,
        ges_fld char(23) not null,
        object char(21) not null,
        filename char(69) not null,
        ges_type char(8) not null,
        setup char(4) not null,
        wg char(4) not null,
        ra real not null,
        dec real not null,
        snr real not null,
        vel real,
        e_vel real,
        vrot real,
        e_vrot real,
        teff real,
        e_teff real,
        nn_teff integer,
        enn_teff real,
        nne_teff integer,
        sys_err_teff real,
        logg real,
        e_logg real,
        nn_logg integer,
        enn_logg real,
        nne_logg integer,
        sys_err_logg real,
        feh real,
        e_feh real,
        nn_feh integer,
        enn_feh real,
        nne_feh integer,
        sys_err_feh real,
        xi real,
        e_xi real,
        nn_xi integer,
        enn_xi real,
        nne_xi integer,
        mh real,
        e_mh real,
        nn_mh integer,
        enn_mh real,
        nne_mh integer,
        alpha_fe real,
        e_alpha_fe real,
        nn_alpha_fe integer,
        enn_alpha_fe real,
        nne_alpha_fe integer,
        vrad real,
        e_vrad real,
        vsini real,
        e_vsini real,
        lim_vsini integer,
        teff_phot real,
        e_teff_phot real,
        teff_irfm_2mj real,
        e_teff_irfm_2mj real,
        teff_irfm_2mh real,
        e_teff_irfm_2mh real,
        teff_irfm_2mks real,
        e_teff_irfm_2mks real,
        teff_irfm_2mave real,
        e_teff_irfm_2mave real,
        fbol_irfm_2m real,
        li1 real,
        upper_combined_li1 integer,
        e_li1 real,
        nn_li1 integer,
        enn_li1 real,
        nl_li1 integer,
        c1 real,
        upper_combined_c1 integer,
        e_c1 real,
        nn_c1 integer,
        enn_c1 real,
        nl_c1 integer,
        c2 real,
        upper_c2 integer,
        e_c2 real,
        nn_c2 integer,
        enn_c2 real,
        nl_c2 integer,
        c3 real,
        upper_c3 integer,
        e_c3 real,
        nn_c3 integer,
        enn_c3 real,
        nl_c3 integer,
        c_c2 real,
        upper_c_c2 integer,
        e_c_c2 real,
        nn_c_c2 integer,
        enn_c_c2 real,
        nl_c_c2 integer,
        n2 real,
        upper_n2 integer,
        e_n2 real,
        nn_n2 integer,
        enn_n2 real,
        nl_n2 integer,
        n3 real,
        upper_n3 integer,
        e_n3 real,
        nn_n3 integer,
        enn_n3 real,
        nl_n3 integer,
        n_cn real,
        upper_n_cn integer,
        e_n_cn real,
        nn_n_cn integer,
        enn_n_cn real,
        nl_n_cn integer,
        o1 real,
        upper_combined_o1 integer,
        e_o1 real,
        nn_o1 integer,
        enn_o1 real,
        nl_o1 integer,
        o2 real,
        upper_o2 integer,
        e_o2 real,
        nn_o2 integer,
        enn_o2 real,
        nl_o2 integer,
        ne1 real,
        upper_combined_ne1 integer,
        e_ne1 real,
        nn_ne1 integer,
        enn_ne1 real,
        nl_ne1 integer,
        ne2 real,
        upper_ne2 integer,
        e_ne2 real,
        nn_ne2 integer,
        enn_ne2 real,
        nl_ne2 integer,
        na1 real,
        upper_combined_na1 integer,
        e_na1 real,
        nn_na1 integer,
        enn_na1 real,
        nl_na1 integer,
        mg1 real,
        upper_combined_mg1 integer,
        e_mg1 real,
        nn_mg1 integer,
        enn_mg1 real,
        nl_mg1 integer,
        al1 real,
        upper_combined_al1 integer,
        e_al1 real,
        nn_al1 integer,
        enn_al1 real,
        nl_al1 integer,
        al3 real,
        upper_al3 integer,
        e_al3 real,
        nn_al3 integer,
        enn_al3 real,
        nl_al3 integer,
        si1 real,
        upper_combined_si1 integer,
        e_si1 real,
        nn_si1 integer,
        enn_si1 real,
        nl_si1 integer,
        si2 real,
        upper_si2 integer,
        e_si2 real,
        nn_si2 integer,
        enn_si2 real,
        nl_si2 integer,
        si3 real,
        upper_si3 integer,
        e_si3 real,
        nn_si3 integer,
        enn_si3 real,
        nl_si3 integer,
        si4 real,
        upper_si4 integer,
        e_si4 real,
        nn_si4 integer,
        enn_si4 real,
        nl_si4 integer,
        s1 real,
        upper_combined_s1 integer,
        e_s1 real,
        nn_s1 integer,
        enn_s1 real,
        nl_s1 integer,
        s2 real,
        upper_s2 integer,
        e_s2 real,
        nn_s2 integer,
        enn_s2 real,
        nl_s2 integer,
        s3 real,
        upper_s3 integer,
        e_s3 real,
        nn_s3 integer,
        enn_s3 real,
        nl_s3 integer,
        ca1 real,
        upper_combined_ca1 integer,
        e_ca1 real,
        nn_ca1 integer,
        enn_ca1 real,
        nl_ca1 integer,
        ca2 real,
        upper_ca2 integer,
        e_ca2 real,
        nn_ca2 integer,
        enn_ca2 real,
        nl_ca2 integer,
        sc1 real,
        upper_combined_sc1 integer,
        e_sc1 real,
        nn_sc1 integer,
        enn_sc1 real,
        nl_sc1 integer,
        sc2 real,
        upper_sc2 integer,
        e_sc2 real,
        nn_sc2 integer,
        enn_sc2 real,
        nl_sc2 integer,
        ti1 real,
        upper_combined_ti1 integer,
        e_ti1 real,
        nn_ti1 integer,
        enn_ti1 real,
        nl_ti1 integer,
        ti2 real,
        upper_ti2 integer,
        e_ti2 real,
        nn_ti2 integer,
        enn_ti2 real,
        nl_ti2 integer,
        v1 real,
        upper_combined_v1 integer,
        e_v1 real,
        nn_v1 integer,
        enn_v1 real,
        nl_v1 integer,
        v2 real,
        upper_v2 integer,
        e_v2 real,
        nn_v2 integer,
        enn_v2 real,
        nl_v2 integer,
        cr1 real,
        upper_combined_cr1 integer,
        e_cr1 real,
        nn_cr1 integer,
        enn_cr1 real,
        nl_cr1 integer,
        cr2 real,
        upper_cr2 integer,
        e_cr2 real,
        nn_cr2 integer,
        enn_cr2 real,
        nl_cr2 integer,
        mn1 real,
        upper_combined_mn1 integer,
        e_mn1 real,
        nn_mn1 integer,
        enn_mn1 real,
        nl_mn1 integer,
        fe1 real,
        upper_combined_fe1 integer,
        e_fe1 real,
        nn_fe1 integer,
        enn_fe1 real,
        nl_fe1 integer,
        fe2 real,
        upper_fe2 integer,
        e_fe2 real,
        nn_fe2 integer,
        enn_fe2 real,
        nl_fe2 integer,
        fe3 real,
        upper_fe3 integer,
        e_fe3 real,
        nn_fe3 integer,
        enn_fe3 real,
        nl_fe3 integer,
        co1 real,
        upper_combined_co1 integer,
        e_co1 real,
        nn_co1 integer,
        enn_co1 real,
        nl_co1 integer,
        ni1 real,
        upper_combined_ni1 integer,
        e_ni1 real,
        nn_ni1 integer,
        enn_ni1 real,
        nl_ni1 integer,
        cu1 real,
        upper_combined_cu1 integer,
        e_cu1 real,
        nn_cu1 integer,
        enn_cu1 real,
        nl_cu1 integer,
        zn1 real,
        upper_combined_zn1 integer,
        e_zn1 real,
        nn_zn1 integer,
        enn_zn1 real,
        nl_zn1 integer,
        sr1 real,
        upper_combined_sr1 integer,
        e_sr1 real,
        nn_sr1 integer,
        enn_sr1 real,
        nl_sr1 integer,
        y1 real,
        upper_combined_y1 integer,
        e_y1 real,
        nn_y1 integer,
        enn_y1 real,
        nl_y1 integer,
        y2 real,
        upper_y2 integer,
        e_y2 real,
        nn_y2 integer,
        enn_y2 real,
        nl_y2 integer,
        zr1 real,
        upper_combined_zr1 integer,
        e_zr1 real,
        nn_zr1 integer,
        enn_zr1 real,
        nl_zr1 integer,
        zr2 real,
        upper_zr2 integer,
        e_zr2 real,
        nn_zr2 integer,
        enn_zr2 real,
        nl_zr2 integer,
        nb1 real,
        upper_combined_nb1 integer,
        e_nb1 real,
        nn_nb1 integer,
        enn_nb1 real,
        nl_nb1 integer,
        mo1 real,
        upper_combined_mo1 integer,
        e_mo1 real,
        nn_mo1 integer,
        enn_mo1 real,
        nl_mo1 integer,
        ru1 real,
        upper_combined_ru1 integer,
        e_ru1 real,
        nn_ru1 integer,
        enn_ru1 real,
        nl_ru1 integer,
        ba2 real,
        upper_ba2 integer,
        e_ba2 real,
        nn_ba2 integer,
        enn_ba2 real,
        nl_ba2 integer,
        la2 real,
        upper_la2 integer,
        e_la2 real,
        nn_la2 integer,
        enn_la2 real,
        nl_la2 integer,
        ce2 real,
        upper_ce2 integer,
        e_ce2 real,
        nn_ce2 integer,
        enn_ce2 real,
        nl_ce2 integer,
        pr2 real,
        upper_pr2 integer,
        e_pr2 real,
        nn_pr2 integer,
        enn_pr2 real,
        nl_pr2 integer,
        nd2 real,
        upper_nd2 integer,
        e_nd2 real,
        nn_nd2 integer,
        enn_nd2 real,
        nl_nd2 integer,
        sm2 real,
        upper_sm2 integer,
        e_sm2 real,
        nn_sm2 integer,
        enn_sm2 real,
        nl_sm2 integer,
        eu2 real,
        upper_eu2 integer,
        e_eu2 real,
        nn_eu2 integer,
        enn_eu2 real,
        nl_eu2 integer,
        gd2 real,
        upper_gd2 integer,
        e_gd2 real,
        nn_gd2 integer,
        enn_gd2 real,
        nl_gd2 integer,
        dy2 real,
        upper_dy2 integer,
        e_dy2 real,
        nn_dy2 integer,
        enn_dy2 real,
        nl_dy2 integer,
        spt char(1),
        veil real,
        e_veil real,
        ew_li real,
        lim_ew_li integer,
        e_ew_li real,
        ewc_li real,
        lim_ewc_li integer,
        e_ewc_li real,
        ew_ha_acc real,
        e_ew_ha_acc real,
        ha10 real,
        e_ha10 real,
        ew_ha_chr real,
        e_ew_ha_chr real,
        fha_chr real,
        e_fha_chr real,
        fwzi real,
        e_fwzi real,
        ew_hb_chr real,
        e_ew_hb_chr real,
        fhb_chr real,
        e_fhb_chr real,
        log_mdot_acc real,
        e_log_mdot_acc real,
        log_l_acc real,
        e_log_l_acc real,
        gamma real,
        e_gamma real,
        convol real,
        e_convol real,
        peculi char(300),
        remark char(300),
        tech char(300)
        );""")

    cursor.execute("ALTER TABLE node_results ADD COLUMN id BIGSERIAL PRIMARY KEY;")
    cursor.close()

    return None


if __name__ == "__main__":


    import psycopg2 as pg
    from astropy.io import fits
    from glob import glob

    connection = pg.connect(dbname="arc")
    create_tables(connection)
    connection.commit()

    from glob import glob
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
    cursor.close()
    
    cursor = connection.cursor()
    from astropy.io import fits
    filenames = glob("data/*.fits")
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
    cursor.execute("""
        create index cname_species_index on line_abundances (cname, element, ion);""")

    cursor.close()


    connection.commit()
    connection.close()

    print("Ingestion complete.")