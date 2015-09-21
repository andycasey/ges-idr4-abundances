#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Flags for line abundances. """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"


import data


def set_bitmask(integers):
    """
    Create a single-value bitmask for the given integers.

    :param integers:
        A list of positive integers.

    :type integers:
        list of int
    """

    if any([0 >= _ for _ in integers]):
        raise ValueError("integers must be positive values")
    
    bitmask = sum([2**i for i in map(int, set(integers))])
    if 0 > bitmask:
        raise ValueError("negative bit mask, likely due to too many high flags")
    return bitmask


def unset_bitmask(bitmask):
    """
    Unpack a single-value bitmask.

    :param bitmask:
        The bitmask value.

    :type bitmask:
        int
    """

    return [i for i in range(1, 1 + int(np.log2(bitmask))) if bitmask >> i & 1]


def retrieve_or_create_line_abundance_flag(database, description):
    """
    Retrieve or create a line abundance flag for the description provided.

    :param database:
        A PostgreSQL database connection.

    :type database:
        :class:`~psycopg2.connection`

    :param description:
        The human-readable description of the flag.

    :type description:
        str

    :returns:
        The given unique id for the created flag.
    """

    result = data.retrieve(database,
        "SELECT id from line_abundance_flags WHERE description = %s",
        (description, ))
    if len(result):
        return result[0][0]

    return create_line_abundance_flag(database, description)


def create_line_abundance_flag(database, description):
    """
    Add a new line abundance flag description to the database.

    :param database:
        A PostgreSQL database connection.

    :type database:
        :class:`~psycopg2.connection`

    :param description:
        The human-readable description of the flag.

    :type description:
        str

    :returns:
        The given unique id for the created flag.
    """

    result = data.retrieve(database, """INSERT INTO line_abundance_flags 
        (description) VALUES (%s) RETURNING id""", (description, ))[0][0]
    database.commit()
    return result


def search_line_abundance_flags(database, keyword):
    """
    Search the descriptions of existing line abundance flags and return a
    dictionary containing the (id, description) pairs.

    :param database:
        A PostgreSQL database connection.

    :type database:
        :class:`~psycopg2.connection`

    :param keyword:
        The keyword to search for in the flag descriptions.

    :type keyword:
        str

    :returns:
        A dictionary containing the flag id as keys, and descriptions as values.
    """

    return dict(data.retrieve(database, """SELECT id, description 
        FROM line_abundance_flags WHERE description::tsvector @@ %s::tsquery""",
        (keyword, )))


def flag_exists(database, flag_id):
    """
    Check whether a flag exists.

    :param database:
        A PostgreSQL database connection.

    :type database:
        :class:`~psycopg2.connection`

    :param id:
        The flag id.

    :type id:
        int

    :returns:
        Whether a flag with the requested id exists or not.
    """

    return data.retrieve(database, """SELECT EXISTS(
        SELECT 1 FROM line_abundance_flags WHERE id = %s)""", (flag_id, ))[0][0]


def update_line_abundance_flag(database, flag_ids, query, values=None):
    """
    Update the flags for some line abundances with the given IDs.

    :param database:
        A PostgreSQL database connection.

    :type database:
        :class:`~psycopg2.connection`

    :param query:
        A SQL query on line_abundances that returns the ids of the relevant lines

    :type query:
        str

    :param flag_ids:
        The unique identifiers of the flags to apply. If an integer is given, it
        must be the bitmask value of all flags to apply. Otherwise, the bitmask
        value will be calculated from the flag ids.

    :type flag_ids:
        int or iterable of ints
    """

    try:
        flag_ids = int(flag_ids)
    except (TypeError, ValueError):
        flag_ids = map(int, flag_ids)
    else:
        flag_ids = [flag_ids]
    
    # Check that all flag IDs actually exist.
    exists = { flag_id: flag_exists(database, flag_id) for flag_id in flag_ids }
    if not all(exists.values()):
        raise ValueError("the following flag ids provided do not exist: {0}"\
            .format(", ".join([str(k) for k, v in exists.items() if not v])))

    # Find the ids that will be affected.
    ids = data.retrieve_table(database, query, values)
    if ids is None: return 0
    
    assert "id" in ids.dtype.names

    # TODO: Update any existing flags.
    values = (set_bitmask(flag_ids), tuple(ids["id"]))
    num_rows =  data.update(database, 
        "UPDATE line_abundances SET flags = %s WHERE id IN %s", values)
    database.commit()
    return num_rows


