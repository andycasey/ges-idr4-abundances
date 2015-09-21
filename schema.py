#!/usr/bin/python

""" Schema management """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"


def create_homogenised_line_abundances_table(database):
    """
    Create a table for the homogenised line abundances.

    :param database:
        A PostgreSQL database connection.

    :type database:
        :class:`~psycopg2.connection`
    """

    # Superfluous column removed:
    #object char(21) not null,
    
    cursor = database.cursor()
    cursor.execute("""CREATE TABLE homogenised_line_abundances (
        cname char(16) not null,
        spectrum_filename_stub char (140) not null,
        wavelength real not null,
        element char(2) not null,
        ion integer not null,
        abundance real,
        e_abundance real,
        upper_abundance int default 0,
        num_measurements int default 0,
        node_bitmask int default 0
        );""")
    cursor.execute("""ALTER TABLE homogenised_line_abundances 
        ADD COLUMN id BIGSERIAL PRIMARY KEY;""")
    cursor.close()

    cursor = database.cursor()
    cursor.execute("""CREATE TABLE homogenised_abundances (
        cname char(16) not null,
        spectrum_filename_stub char(140) not null,
        element char(2) not null,
        ion integer not null,
        abundance real,
        e_abundance real,
        upper_abundance int default 0,
        num_lines int default 0,
        num_measurements int default 0);""")
    cursor.execute("""ALTER TABLE homogenised_abundances ADD COLUMN id
        BIGSERIAL PRIMARY KEY;""")
    cursor.close()

    database.commit()
    return True


def drop_homogenised_line_abundances_table(database):
    """
    Drop the table used for homogenised line abundances.

    :param database:
        A PostgreSQL database connection.

    :type database:
        :class:`~psycopg2.connection`
    """

    cursor = database.cursor()
    cursor.execute("DROP TABLE IF EXISTS homogenised_line_abundances;")
    cursor.execute("DROP TABLE IF EXISTS homogenised_abundances;")
    cursor.close()
    database.commit()

    return True