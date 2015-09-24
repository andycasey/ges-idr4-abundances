
import os
import logging
from glob import glob
from collections import OrderedDict
import matplotlib.pyplot as plt

import release

logger = logging.getLogger("ges")


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


def make_figures(figures, element, ion, format="png"):

    directory = "figures/{0}{1}".format(element.upper(), ion)
    if not os.path.exists(directory): os.mkdir(directory)
    for figure_name, details in figures.items():
        
        if not isinstance(details, (tuple, )):
            details = (details, )

        command = details[0]
        kwargs = {} if len(details) == 1 else details[1]
        

        print("Doing {0} with {1}".format(figure_name, kwargs))
        try:
            fig = command(element, ion, **kwargs)
        except:
            logger.exception("Something happened")
            raise
            
        if fig is not None:
        
            if isinstance(fig, dict):
                for suffix, figure in fig.items():
                    filename = os.path.join(directory, "{0}-{1}.{2}".format(
                        figure_name, suffix, format))
                    if figure is not None:
                        savefig(figure, filename)
                        print("Created {}".format(filename))

            else:
                filename = os.path.join(directory, "{0}.{1}".format(
                    figure_name, format))
                savefig(fig, filename)
                print("Created {}".format(filename))


    # Create an all.html
    filenames = glob("{0}/*.{1}".format(directory, format))
    html = "<html><body>"
    for filename in filenames:
        html += "{0}: <a href=\"{0}\"><img src=\"{0}\" /></a><br />".format(
            os.path.basename(filename))

    html += "</body></html>"
    with open("{0}/all.html".format(directory), "w") as fp:
        fp.write(html)



ges = release.DataRelease("arc")

"""

"""
species = [
    ("Si", 1, None),
    ("Si", 2, None),
    ("Al", 1, (3, 7.5)),
    ("Al", 3, (3, 7.5)),
    ("Na", 1, (2.5, 7.5)),
    ("S", 1, None),
    ("S", 2, None),
    ("S", 3, None),
    ("Ne", 1, None),
    ("Ne", 2, None),
    
    ("Mg", 1, (5, 9)),
    ("Ca", 1, (3.5, 8)),
    ("Ca", 2, (3.5, 8)),

    ("Ni", 1, (3, 10)),
    ("Zn", 1, (2, 7)),
    ("Cr", 1, (0, 10)),
    ("Cr", 2, (0, 10)),
    ("Sc", 1, (0, 6)),
    ("Sc", 2, (0, 6)),
    ("Cu", 1, (1, 6)),
    ("Co", 1, (0, 6)),
    ("Mn", 1, None),
    ("Zr", 1, None),
    ("Zr", 2, None),
    ("V", 1, None),
    ("V", 2, None),
    ("O", 1, None),
    ("O", 2, None),
    ("Ti", 1, None),
    ("Ti", 2, None),
    ("Fe", 1, None),
    ("Fe", 2, None),
    ("Fe", 3, None),
    ("La", 2, None),
    ("Ba", 2, None),
    ("Y", 1, None),
    ("Y", 2, None),
    ("C", 1, None),
    ("C", 2, None),
    ("C", 3, None),
    #("C_C", )
    #("N_C", )
    ("N", 1, None),
    ("N", 2, None),
    ("Li", 1, None),
    ("Ru", 1, None),
    ("Sr", 1, None),
    ("Dy", 2, None),
    ("Ce", 2, None),
    ("Pr", 2, None),
    ("Nb", 1, None),
    ("Nd", 2, None),
    ("Mo", 1, None),
    ("Eu", 2, None),
    ("Gd", 2, None),
    ("Sm", 2, None)
]

for element, ion, absolute_extent in species:
    figures = OrderedDict([
        ("differential-line-abundances-wrt-teff", (
            ges.plot.differential_line_abundances_wrt_parameter,
            { "parameter": "teff", "x_extent": (3500, 7000) })),


    ])

    figures = OrderedDict([

        ("differential-line-abundances-wrt-teff", (
            ges.plot.differential_line_abundances_wrt_parameter,
            { "parameter": "teff", "x_extent": (3500, 7000) })),
        ("differential-line-abundances-wrt-logg", (
            ges.plot.differential_line_abundances_wrt_parameter,
            { "parameter": "logg", "x_extent": (0, 5) })),
        ("differential-line-abundances-wrt-feh", (
            ges.plot.differential_line_abundances_wrt_parameter,
            { "parameter": "feh", "x_extent": (-3, 0.5) })),

        ("compare-bm", (ges.plot.benchmark_comparison,
            { "benchmark_filename": "benchmarks.yaml" })),
        ("compare-solar", ges.plot.solar_comparison),
        ("compare-m67-1194", ges.plot.m67_twin_comparison),
        ("benchmarks", (ges.plot.benchmark_line_abundances, { 
            "benchmark_filename": "benchmarks.yaml" })),
        ("differential-line-abundances", ges.plot.differential_line_abundances),
        ("differential-line-abundances-clipped", (
            ges.plot.differential_line_abundances, { "absolute_extent": absolute_extent })),
        ("abundance-heatmap", (ges.plot.transition_heatmap, {"column": "abundance"})),
        ("ew-heatmap", (ges.plot.transition_heatmap, {"column": "ew"})),
        ("mean-abundances", (ges.plot.mean_abundances, {
            "extent": absolute_extent })),

        ("line-abundances-logx-wrt-teff", (ges.plot.line_abundances, {
            "reference_column": "teff",
            "abundance_format": "log_x",
            "aux_column": "logg",
            "aux_extent": (0, 5),
            "extent": None,
            "highlight_flagged": True,
            "show_node_comparison": True,
            "show_line_comparison": True,
            "show_homogenised": True,
            })),
        ("line-abundances-logx-wrt-logg", (ges.plot.line_abundances, {
            "abundance_format": "log_x",
            "reference_column": "logg",
            "aux_column": "feh",
            "aux_extent": (-3, 1),
            "extent": None,
            "highlight_flagged": True,
            "show_node_comparison": True,
            "show_line_comparison": True,
            "show_homogenised": True,
            })),
        ("line-abundances-logx-wrt-feh", (ges.plot.line_abundances, {
            "abundance_format": "log_x",
            "reference_column": "feh",
            "aux_column": "teff",
            "aux_extent": (3500, 7000),
            "extent": None,
            "highlight_flagged": True,
            "show_node_comparison": True,
            "show_line_comparison": True,
            "show_homogenised": True,
            })),
        ("line-abundances-logx-wrt-teff-clip", (ges.plot.line_abundances, {
            "abundance_format": "log_x",
            "reference_column": "teff",
            "aux_column": "logg",
            "aux_extent": (0, 5),
            "extent": absolute_extent,
            "highlight_flagged": True,
            "show_node_comparison": True,
            "show_line_comparison": True,
            "show_homogenised": True,
            })),
        ("line-abundances-logx-wrt-logg-clip", (ges.plot.line_abundances, {
            "abundance_format": "log_x",
            "reference_column": "logg",
            "aux_column": "feh",
            "aux_extent": (-3, 1),
            "extent": absolute_extent,
            "highlight_flagged": True,
            "show_node_comparison": True,
            "show_line_comparison": True,
            "show_homogenised": True,
            })),
        ("line-abundances-logx-wrt-feh-clip", (ges.plot.line_abundances, {
            "abundance_format": "log_x",
            "reference_column": "feh",
            "aux_column": "teff",
            "aux_extent": (3500, 7000),
            "extent": absolute_extent,
            "highlight_flagged": True,
            "show_node_comparison": True,
            "show_line_comparison": True,
            "show_homogenised": True,
            })),
        ("line-abundances-xfe-wrt-teff", (ges.plot.line_abundances, {
            "abundance_format": "x_fe",
            "reference_column": "teff", "aux_column": "logg",
            "aux_extent": (0, 5) })),
        ("line-abundances-xfe-wrt-logg", (ges.plot.line_abundances, {
            "abundance_format": "x_fe",
            "reference_column": "logg", "aux_column": "feh",
            "aux_extent": (-3, 1) })),
        ("line-abundances-xfe-wrt-feh", (ges.plot.line_abundances, {
            "abundance_format": "x_fe",
            "reference_column": "feh", "aux_column": "teff",
            "aux_extent": (3500, 7000) })),


    ])

    
    try:
        make_figures(figures, element, ion)
    except:
        logger.exception("FAIL SNAIL")
        raise 

    plt.close("all")