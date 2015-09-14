#!/usr/bin/python

""" Plots."""

import logging
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from colormaps import (magma, inferno, plasma, viridis)

import numpy as np
from astropy.table import Table



def label(text_label):
    return text_label.replace("_", "-")


def _scatter_elements(ax, node, x_label, y_label, limits=True, errors=True,
    color_label=None, **kwargs):

    # TODO Show upper limits correctly.
    # TODO show those with flags correctly
    # TODO show uncertainties.


    elements = [column for _, column, comment in node[2].header.cards \
        if "Abundance" in comment]

    kwds = {
        "c": "#666666",
        "cmap": plasma
    }
    kwds.update(kwargs)
    if color_label is not None:
        kwds["c"] = node[2].data[color_label]

    return ax.scatter(node[2].data[x_label], node[2].data[y_label], **kwds)




def node_element_abundance_against_stellar_parameters(node, element,
    limits=True, errors=True):
    """
    Plot the reported node abundance against stellar parameters.

    :param node:
        The node data.

    :type node:
        FITS image

    :param element:
        The name of the element to show.

    :type element:
        str
    """


    element = element.upper()

    data = node[2].data

    fig, (ax_teff, ax_logg, ax_feh) = plt.subplots(3)


    _scatter_elements(ax_teff, node, "TEFF", element,
        limits=limits, errors=errors, color_label="NL_{}".format(element))

    _scatter_elements(ax_logg, node, "LOGG", element,
        limits=limits, errors=errors, color_label="NL_{}".format(element))

    scat = _scatter_elements(ax_feh, node, "FEH", element,
            limits=limits, errors=errors, color_label="NL_{}".format(element))

    # Labels
    ax_teff.set_xlabel(label("TEFF"))
    ax_logg.set_xlabel(label("LOGG"))
    ax_feh.set_xlabel(label("FEH"))
    [ax.set_ylabel(label(element)) for ax in fig.axes]

    for ax in fig.axes:
        ax.xaxis.set_major_locator(MaxNLocator(6))
        ax.yaxis.set_major_locator(MaxNLocator(6))

    # Title.
    ax_teff.set_title("GES {release} {node} ({instrument} {date})".format(
        release=node[0].header.get("RELEASE", "iDR4?"),
        date=node[0].header.get("DATETAB", "UNKNOWN DATE"),
        instrument=node[0].header.get("INSTRUME", "UNKNOWN INSTRUMENT"),
        node=node[0].header["NODE1"]))

    fig.tight_layout()

    # Colorbar
    ax_cbar = fig.colorbar(scat, ax=fig.axes)
    ax_cbar.set_label(label("NL_{}".format(element)))

    return fig


def retrieve_table(database, query, values=None):

    cursor = database.cursor()
    cursor.execute(query, values)
    names = [_[0] for _ in cursor.description]
    results = cursor.fetchall()
    if len(results) > 0:
        t = Table(rows=results, names=names)
        cursor.close()
        return t
    return None


def transition_heatmap(database, element, ion, column="abundance", linear=False,
    **kwargs): 
    """
    Display a heatmap of lines that were used for a given species.

    :param database:
        A PostgreSQL database connection.

    :param element:
        The name of the element to display.

    :type element:
        str

    :param ion:
        The ionisation state of the element (1 = neutral).

    :type ion:
        int

    :param column: [optional]
        The column name to display the heat map for. Default is abundance.

    :type column:
        str
    """

    # Retrieve all line abundances for this element and ion.
    data = retrieve_table(database,
        "SELECT * FROM line_abundances WHERE element = %s and ion = %s",
        (element, ion))
    if data is None: return

    column = column.lower()
    if column not in data.dtype.names:
        raise ValueError("column '{0}' does not exist".format(column))

    # Identify the unique nodes and wavelengths.
    nodes = sorted(set(data["node"]))
    wavelengths = sorted(set(data["wavelength"]))
    N_nodes, N_wavelengths = map(len, (nodes, wavelengths))

    # Build a count map array.
    count = np.zeros((N_nodes, N_wavelengths))
    for i, node in enumerate(nodes):
        for j, wavelength in enumerate(wavelengths):
            count[i, j] = np.sum(np.isfinite(data[column][
                (data["wavelength"] == wavelength) * (data["node"] == node) \
                    * (data["upper_{}".format(column)] == 0)]))
    kwds = {
        "aspect": "auto",
        "cmap": plasma,
        "interpolation": "nearest"
    }
    kwds.update(kwargs)

    fig, ax = plt.subplots(
        figsize=(6.5 + N_wavelengths * 0.25, 2 + N_nodes * 0.25))
    if linear:
        # Create the wavelength map.
        px = np.diff(wavelengths).min()
        wavelength_map = np.arange(min(wavelengths), max(wavelengths) + px, px)

        # Fill in the values.
        heat_map = np.nan * np.ones((N_nodes, len(wavelength_map), 1))
        for i in range(N_nodes):
            for j, wavelength in enumerate(wavelengths):
                index = wavelength_map.searchsorted(wavelength)
                heat_map[i, index, :] = count[i, j]


        raise NotImplementedError

    else:
        image = ax.imshow(count, **kwds)

    ax.set_xlabel(r"Wavelength $[\AA]$")
    ax.set_title("{element} {ion} (measured {column})".format(element=element,
        ion=ion, column=column))

    ax.set_xticks(np.arange(N_wavelengths) - 0.5)
    ax.set_xticklabels(["{0:.1f}".format(_) for _ in wavelengths],
        rotation=45)

    ax.set_yticks(np.arange(N_nodes))
    ax.set_yticklabels(nodes)

    ax.xaxis.set_tick_params(width=0)
    ax.yaxis.set_tick_params(width=0)


    fig.tight_layout()

    cbar = plt.colorbar(image, ax=[ax])
    cbar.locator = MaxNLocator(3)
    cbar.update_ticks()
    cbar.ax.set_aspect(2)
    cbar.set_label(r"$N$")


    return fig


def tellurics(database, vel_bin=1, wavelength_bin=0.5,
    vel_range=(-300, 300)):

    #line_data = retrieve_table(line_database,
    #    "SELECT * FROM line_abundances WHERE element = %s and ion = %s",
    #    (element, ion))
    line_data = retrieve_table(database, "SELECT * FROM line_abundances")
    if line_data is None: return
    N_lines = len(line_data)

    # Create a map for the velocity and wavelength bins.
    wavelengths = np.arange(line_data["wavelength"].min(),
        line_data["wavelength"].max() + wavelength_bin, wavelength_bin)

    # Get all the velocities.
    # Just get one of the nodes.
    cursor = database.cursor()
    cursor.execute("SELECT DISTINCT(node) from node_results LIMIT 1")
    node = cursor.fetchone()
    cursor.close()

    vel_data = retrieve_table(database,
        "SELECT cname, vel, node FROM node_results WHERE node = %s", (node, ))

    vel_range = [None, None] if vel_range is None else list(vel_range)
    if vel_range[0] is None: vel_range[0] = np.nanmin(vel_data["vel"])
    if vel_range[1] is None: vel_range[1] = np.nanmax(vel_data["vel"])

    velocities = np.arange(vel_range[0], vel_range[1], vel_bin)

    abundance_stddev = np.zeros((wavelengths.size, velocities.size))
    N = abundance_stddev.copy()

    # For each abundance line, match to CNAME, find which bin it is in,
    # calculate.

    for i, line in enumerate(line_data):

        print("{0}/{1}".format(i, N_lines))

        # Match to cname
        v_index = np.where(vel_data["cname"] == line["cname"])[0]
        assert len(v_index) == 1

        i_index = wavelengths.searchsorted(line["wavelength"]) - 1
        j_index = velocities.searchsorted(vel_data["vel"][v_index[0]]) - 1

        # OK, now that we know where it is, let's calculate the abundance offset
        # for this line with respect to all other lines of the same species 
        comparison_data = retrieve_table(database, 
            """SELECT abundance FROM line_abundances WHERE element = %s and 
            ion = %s and filename = %s and upper_abundance = 0""",
            (line["element"], line["ion"], line["filename"]))

        if len(comparison_data) == 1: continue        

        mean = np.nanmean(comparison_data["abundance"])
        sigma = np.nanstd(comparison_data["abundance"])

        abundance_stddev[i_index, j_index] += (line["abundance"] - mean)/sigma
        N[i_index, j_index] += 1



    raise a


def abundance_covariance(database, element, ion, node=None):
    """
    Show the covariance in all line abundances for a given species.
    """

    if node is None:
        query_suffix, args_suffix = "", []
    else:
        query_suffix, args_suffix = " AND node = %s", [node]


    args = [element, ion] + args_suffix
    data = retrieve_table(database, """SELECT wavelength, filename, abundance 
        FROM line_abundances WHERE element = %s and ion = %s""" + query_suffix,
         args)
    if data is None: return None

    # Arrange everything by filename.
    filenames = retrieve_table(database, """SELECT DISTINCT(filename) FROM 
        line_abundances WHERE element = %s AND ion = %s""" + query_suffix, args)
    wavelengths = retrieve_table(database, """SELECT DISTINCT(wavelength) FROM
        line_abundances WHERE element = %s AND ion = %s""" + query_suffix, args)
    filenames, wavelengths = filenames["filename"], wavelengths["wavelength"]

    extra_matches = {}
    std_devs = {}
    X = np.nan * np.ones((len(filenames), len(wavelengths)))
    for i, filename in enumerate(filenames):
        for j, wavelength in enumerate(wavelengths):
            #print(i * len(wavelengths) + j, X.size)

            indices = (data["filename"] == filename) \
                * (data["wavelength"] == wavelength)
            if indices.sum() > 0:
                if indices.sum() > 1:
                    print("Warning: {0} matches: {1}".format(indices.sum(),
                        data["abundance"][indices]))
                    _ = "{0}.{1}".format(data["filename"], data["wavelength"])
                    extra_matches[_] = indices.sum()
                    std_devs[_] = np.nanstd(data["abundance"][indices])
                X[i, j] = np.nanmean(data["abundance"][indices])


    raise a
    #line_data = retrieve_table(line_database,
    #    "SELECT * FROM line_abundances WHERE element = %s and ion = %s",
    #    (element, ion))


if __name__ == "__main__":


    """
    from data import load_node

    d = load_node("data/GES_iDR4_WG11_MyGIsFOS.fits")
    f = node_element_abundance_against_stellar_parameters(d, "SI1")
    """

    import psycopg2 as pg
    db = pg.connect(dbname="arc")

    #tellurics(db)#, "Si", 1)
    abundance_covariance(db, "Si", 1)

    transition_heatmap(db, "Si", 1)
    transition_heatmap(db, "Si", 2)
    #transition_heatmap(db, "Fe", 1)
    #transition_heatmap(db, "Fe", 2)
    raise a