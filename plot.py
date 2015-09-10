#!/usr/bin/python

""" Plots."""

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from colormaps import (magma, inferno, plasma, viridis)


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
        release=node[0].header["RELEASE"], date=node[0].header["DATETAB"],
        instrument=node[0].header["INSTRUME"], node=node[0].header["NODE1"]))

    fig.tight_layout()

    # Colorbar
    ax_cbar = fig.colorbar(scat, ax=fig.axes)
    ax_cbar.set_label(label("NL_{}".format(element)))

    return fig



if __name__ == "__main__":

    from data import load_node

    d = load_node("data/GES_iDR4_WG11_MyGIsFOS.fits")

    f = node_element_abundance_against_stellar_parameters(d, "SI1")

    raise a