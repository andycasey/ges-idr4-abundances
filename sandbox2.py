
import os
import logging
from glob import glob
from collections import OrderedDict
import matplotlib.pyplot as plt

import release

logger = logging.getLogger("ges")

wg = "WG10"


def savefig(figure, filename):

    dpi = plt.rcParams["savefig.dpi"]

    while True:
        try:
            figure.savefig(filename, dpi=dpi)
        except ValueError:
            dpi /= 2.
            logger.warn("Dropping DPI to {0} to save figure {1}".format(dpi,
                filename))
        else:
            break

    return True


def make_figures(figures, element, ion, format="png", remake=True):

    failed = []
    directory = "figures/{0}/{1}{2}".format(wg, element.upper(), ion)
    if not os.path.exists(directory): os.mkdir(directory)
    for figure_name, details in figures.items():
        
        if not isinstance(details, (tuple, )):
            details = (details, )

        command = details[0]
        kwargs = {} if len(details) == 1 else details[1]



        expected_filename = "{0}/{1}.{2}".format(directory, figure_name, format)

        if not remake and os.path.exists(expected_filename):
            logger.info("Figure already exists: {0},.. skipping".format(
                expected_filename))
            continue
    

        print("Doing {0} with {1}".format(figure_name, kwargs))
        try:
            fig = command(element, ion, **kwargs)
        except:
            logger.exception("Something happened")
            failed.append((figure_name, element, ion))
            raise
            continue

        if fig is not None:
        
            if isinstance(fig, dict):
                for suffix, figure in fig.items():
                    filename = os.path.join(directory, "{0}-{1}.{2}".format(
                        figure_name, suffix, format))
                    if figure is not None:
                        savefig(figure, filename)
                        print("Created {}".format(filename))
                        plt.close(figure)
            else:
                filename = os.path.join(directory, "{0}.{1}".format(
                    figure_name, format))
                savefig(fig, filename)
                print("Created {}".format(filename))
                plt.close(fig)

            plt.close("all")


    # Create an all.html
    filenames = glob("{0}/*.{1}".format(directory, format))
    html = "<html><body>"
    for filename in filenames:
        html += "{0}: <a href=\"{0}\"><img src=\"{0}\" /></a><br />".format(
            os.path.basename(filename))

    html += "</body></html>"
    with open("{0}/all.html".format(directory), "w") as fp:
        fp.write(html)

    return failed


ges = release.DataRelease("ges-idr4-wg10")

species = [ # As taken from the FITS table headers.
    ("Li", 1, None),
    ("C", 1, None),
    ("C", 2, None),
    ("C", 3, None),
    ("C_C2", 1, None),
    ("N", 2, None),
    ("N", 3, None),
    ("N_CN", 1, None),
    ("O", 1, None),
    ("O", 2, None),
    ("Ne", 1, None),
    ("Ne", 2, None),
    ("Na", 1, None),
    ("Mg", 1, None),
    ("Al", 1, None),
    ("Al", 3, None),
    ("Si", 1, None),
    ("Si", 2, None),
    ("Si", 3, None),
    ("Si", 4, None),
    ("S", 1, None),
    ("S", 2, None),
    ("S", 3, None),
    ("Ca", 1, None),
    ("Ca", 2, None),
    ("Sc", 1, None),
    ("Sc", 2, None),
    ("Ti", 1, None),
    ("Ti", 2, None),
    ("V", 1, None),
    ("V", 2, None),
    ("Cr", 1, None),
    ("Cr", 2, None),
    ("Mn", 1, None),
    ("Fe", 1, None),
    ("Fe", 2, None),
    ("Fe", 3, None),
    ("Co", 1, None),
    ("Ni", 1, None),
    ("Cu", 1, None),
    ("Zn", 1, None),
    ("Sr", 1, None),
    ("Y", 1, None),
    ("Y", 2, None),
    ("Zr", 1, None),
    ("Zr", 2, None),
    ("Nb", 1, None),
    ("Mo", 1, None),
    ("Ru", 1, None),
    ("Ba", 2, None),
    ("La", 2, None),
    ("Ce", 2, None),
    ("Pr", 2, None),
    ("Nd", 2, None),
    ("Sm", 2, None),
    ("Eu", 2, None),
    ("Gd", 2, None),
    ("Dy", 2, None),
]



for element, ion, absolute_extent in species:
    

    figures = OrderedDict([
        ("mean-abundances", (ges.plot.mean_abundances, {
    ])

    try:
        make_figures(figures, element, ion)
    except:
        logger.exception("FAIL SNAIL")
        raise 

    plt.close("all")

