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


def _scatter_elements(ax, data, x_label, y_label, limits=True, errors=True,
    color_label=None, **kwargs):

    # TODO Show upper limits correctly.
    # TODO show those with flags correctly
    # TODO show uncertainties.


    kwds = {
        "c": "#666666",
        "cmap": plasma
    }
    kwds.update(kwargs)
    if color_label is not None:
        kwds["c"] = data[color_label]

    return ax.scatter(data[x_label], data[y_label], **kwds)


def mean_abundance_against_stellar_parameters(database, element, ion, node=None,
    limits=True, errors=True, **kwargs):
    """
    Plot the reported node abundance against stellar parameters.

    """

    if node is None:
        # Get all the node names.
        nodes = list(retrieve_table(database,
            "SELECT DISTINCT(node) from node_results")["node"])
    else:
        nodes = [node]


    figures = {}    
    column = "{0}{1}".format(element.lower(), ion)

    for node in nodes:

        fig, (ax_teff, ax_logg, ax_feh) = plt.subplots(3)
        data = retrieve_table(database, "SELECT teff, logg, feh, {0}, nl_{0} "\
            "FROM node_results WHERE node = %s".format(column), (node, ))

        _scatter_elements(ax_teff, data, "teff", column, limits=limits,
            errors=errors, color_label="nl_{}".format(column))

        _scatter_elements(ax_logg, data, "logg", column, limits=limits,
            errors=errors, color_label="nl_{}".format(column))
        
        scat = _scatter_elements(ax_feh, data, "feh", column, limits=limits,
            errors=errors, color_label="nl_{}".format(column))
        
        # Labels
        ax_teff.set_xlabel(label("TEFF"))
        ax_logg.set_xlabel(label("LOGG"))
        ax_feh.set_xlabel(label("FEH"))
        [ax.set_ylabel(label(column.upper())) for ax in fig.axes]

        for ax in fig.axes:
            ax.xaxis.set_major_locator(MaxNLocator(6))
            ax.yaxis.set_major_locator(MaxNLocator(6))

        # Title.
        wg = retrieve_table(database, "SELECT wg from node_results limit 1")
        ax_teff.set_title("GES {release} {wg} {node}".format(
            release=kwargs.pop("release", "iDR4"), node=node.strip(),
            wg=wg["wg"][0]))

        fig.tight_layout()

        # Colorbar
        ax_cbar = fig.colorbar(scat, ax=fig.axes)
        ax_cbar.set_label(label("NL_{}".format(column.upper())))

        figures[node.strip()] = fig

    return figures


def retrieve(database, query, values=None):
    cursor = database.cursor()
    cursor.execute(query, values)
    result = cursor.fetchall()
    cursor.close()
    return result


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
        #assert len(v_index) == 1

        i_index = wavelengths.searchsorted(line["wavelength"]) - 1
        j_index = velocities.searchsorted(np.nanmean(vel_data["vel"][v_index])) - 1

        # OK, now that we know where it is, let's calculate the abundance offset
        # for this line with respect to all other lines of the same species 

        indices = \
            (line_data["filename"] == line["filename"]) * \
            (line_data["element"] == line["element"]) * \
            (line_data["ion"] == line["ion"]) * \
            (line_data["upper_abundance"] == 0)

        comparison_data = line_data[indices]


        #comparison_data = retrieve_table(database, 
        #    """SELECT abundance FROM line_abundances WHERE element = %s and 
        #    ion = %s and filename = %s and upper_abundance = 0""",
        #    (line["element"], line["ion"], line["filename"]))

        if len(comparison_data) == 1: continue        

        mean = np.nanmean(comparison_data["abundance"])
        sigma = np.nanstd(comparison_data["abundance"])

        abundance_stddev[i_index, j_index] += (line["abundance"] - mean)/sigma
        N[i_index, j_index] += 1



    raise a


def transition_covariance(database, element, ion, node=None, column="abundance",
    vmin=None, vmax=None, **kwargs):
    """
    Show the covariance in all line abundances for a given species.
    """

    if node is None:
        query_suffix, args_suffix = "", []
    else:
        query_suffix, args_suffix = " AND node = %s", [node]

    column = column.lower()
    if column not in ("ew", "abundance"):
        raise ValueError("column must be ew or abundance")

    args = [element, ion] + args_suffix
    data = retrieve_table(database, """SELECT wavelength, filename, ew, abundance 
        FROM line_abundances WHERE element = %s and ion = %s""" + query_suffix,
         args)
    if data is None: return None

    # Arrange everything by filename.
    filenames = retrieve_table(database, """SELECT DISTINCT(filename) FROM 
        line_abundances WHERE element = %s AND ion = %s""" + query_suffix, args)
    wavelengths = retrieve_table(database, """SELECT DISTINCT(wavelength) FROM
        line_abundances WHERE element = %s AND ion = %s""" + query_suffix, args)
    
    filenames = filenames["filename"]
    wavelengths = np.sort(wavelengths["wavelength"])

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
                        data[column][indices]))
                    _ = "{0}.{1}".format(filename, wavelength)
                    extra_matches[_] = indices.sum()
                    std_devs[_] = np.nanstd(data[column][indices])
                X[i, j] = np.nanmean(data[column][indices])

    # If there are *no* finite measurements in a row, skip it.
    X = X[np.any(np.isfinite(X), axis=1)]
    X = np.ma.array(X, mask=~np.isfinite(X))
    cov = np.ma.cov(X.T)

    kwds = {
        "aspect": "auto",
        "cmap": plasma,
        "interpolation": "nearest",
        "vmin": vmin,
        "vmax": vmax
    }
    kwds.update(kwargs)
    fig, ax = plt.subplots()

    ax.patch.set_facecolor("#CCCCCC")
    image = ax.imshow(cov, **kwds)
    
    _ = np.arange(len(wavelengths)) + 0.5
    ax.set_xticks(_, minor=True)
    ax.set_yticks(_, minor=True)
    ax.set_xticks(_)
    ax.set_yticks(_ + 0.5)
    ax.xaxis.set_tick_params(width=0, which="major")
    ax.yaxis.set_tick_params(width=0, which="major")
    ax.set_xticklabels(["{0:.1f}".format(_) for _ in wavelengths], rotation=45)
    ax.set_yticklabels(["{0:.1f}".format(_) for _ in wavelengths], rotation=0)
    
    ax.set_xlim(0.5, len(wavelengths) + 0.5)
    ax.set_ylim(0.5, len(wavelengths) + 0.5)
    ax.set_xlabel(r"Wavelength $[\AA]$")
    ax.set_ylabel(r"Wavelength $[\AA]$")

    ax.set_title("{node} {element} {ion} (measured {column})".format(
        element=element, ion=ion, column=column, node=node or "All nodes"))

    fig.tight_layout()
    cbar = plt.colorbar(image, ax=[ax])
    return fig





def corner_scatter(data, labels=None, uncertainties=None, extent=None,
    color=None, bins=20):
    """
    Create a corner scatter plot showing the differences between each node.

    :param extent: [optional]
        The (minimum, maximum) extent of the plots.

    :type extent:
        two-length tuple

    :returns:
        A matplotlib figure.
    """

    # How many nodes to plot?
    N = data.shape[0]
    K = N
    assert K > 0, "Need more than one node to compare against."

    factor = 2.0           # size of one side of one panel
    lbdim = 0.5 * factor   # size of left/bottom margin
    trdim = 0.5 * factor   # size of top/right margin
    whspace = 0.15         # w/hspace size
    plotdim = factor * K + factor * (K - 1.) * whspace
    dim = lbdim + plotdim + trdim

    fig, axes = plt.subplots(K, K, figsize=(dim, dim))
    
    lb = lbdim / dim
    tr = (lbdim + plotdim) / dim
    fig.subplots_adjust(left=lb, bottom=lb, right=tr, top=tr,
        wspace=whspace, hspace=whspace)

    hist_kwargs = {
        "color": "k",
        "histtype": "step"
    }

    if extent is None:
        extent = (0.9 * np.nanmin(data), 1.1 * np.nanmax(data))
    
    # Match all of the nodes
    for i in range(N):
        for j in range(N):

            if j > i: # hide.
                try:
                    ax = axes[i, j]
                    ax.set_visible(False)
                    ax.set_frame_on(False)
                except IndexError:
                    None
                continue
                
            ax = axes[i, j]

            if i == j:
                indices = np.arange(N)
                indices = indices[indices != i]

                diff = (data[i] - data[indices]).flatten()
                diff = diff[np.isfinite(diff)]
                if diff.size:
                    ax.hist(diff, bins=bins, **hist_kwargs)
                
            else:    
                ax.plot(extent, extent, "k:", zorder=-100)
                if uncertainties is not None:
                    ax.errorbar(data[i], data[j], 
                        xerr=uncertainties[i], yerr=uncertainties[j],
                        ecolor="k", aa=True, fmt=None, mec='k', mfc="w", ms=6,
                        zorder=1, lc="k")

                ax.scatter(data[i], data[j], c=color or "w", zorder=100)
                
                ax.set_xlim(extent)
                ax.set_ylim(extent)

            ax.xaxis.set_major_locator(MaxNLocator(5))
            ax.yaxis.set_major_locator(MaxNLocator(5))

            if i < K - 1:
                ax.set_xticklabels([])
            else:
                [l.set_rotation(45) for l in ax.get_xticklabels()]
                if labels is not None:
                    if i == j:
                        ax.set_xlabel("{} $-$ $X$".format(labels[j]))
                    else:
                        ax.set_xlabel(labels[j])
                    ax.xaxis.set_label_coords(0.5, -0.3)
                
            if j > 0 or i == j:
                ax.set_yticklabels([])
            else:
                [l.set_rotation(45) for l in ax.get_yticklabels()]
                if labels is not None:
                    ax.set_ylabel(labels[i])
                    ax.yaxis.set_label_coords(-0.3, 0.5)


    #fig.tight_layout()

    return fig


def mean_abundance_differences(database, element, ion, bins=None):
    """
    Show the mean abundance differences from each node.
    """

    nodes = retrieve_table(database, 
        "SELECT distinct(node) from node_results")["node"]
    N_entries = int(retrieve(database, """SELECT count(*) FROM node_results
        WHERE node = %s""", (nodes[0], ))[0][0])
    N_nodes = len(nodes)

    Z = np.nan * np.ones((N_nodes, N_entries))
    Z_uncertainties = Z.copy()
    #data_upper = data.copy()
    column = "{0}{1}".format(element.lower(), ion)
    for i, node in enumerate(nodes):

        data = retrieve_table(database,
            """SELECT {0}, e_{0} FROM node_results WHERE node = %s
            ORDER BY ra DESC, dec DESC, cname DESC""".format(column), (node, ))
        Z[i, :] = data[column]
        Z_uncertainties[i, :] = data["e_{}".format(column)]
        #data_upper[i, :] = data["upper_{}".format(column)]

    bins = bins or np.arange(-0.50, 0.55, 0.05)
    return corner_scatter(Z, uncertainties=Z_uncertainties,
        labels=map(str.strip, nodes), bins=bins)


    # x-axis = average REW
    # y-axis = node1 - node2




if __name__ == "__main__":


    """
    from data import load_node

    d = load_node("data/GES_iDR4_WG11_MyGIsFOS.fits")
    f = mean_abundances_against_stellar_parameters(d, "SI1")
    """

    import psycopg2 as pg
    db = pg.connect(dbname="arc")

    #tellurics(db)#, "Si", 1)
    #transition_covariance(db, "Si", 1)
    #mean_abundance_against_stellar_parameters(db, "Si", 1)
    mean_abundance_differences(db, "Si", 1)
    raise a
    #

    transition_heatmap(db, "Si", 1)
    transition_heatmap(db, "Si", 2)
    #transition_heatmap(db, "Fe", 1)
    #transition_heatmap(db, "Fe", 2)
    raise a