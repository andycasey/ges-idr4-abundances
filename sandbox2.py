
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


def make_figures(figures, element, ion, format="png", remake=True):

    failed = []
    directory = "figures/{0}{1}".format(element.upper(), ion)
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


ges = release.DataRelease("arc")

species = [
    #("Al", 1, (3, 7.5)),
    #("Ba", 2, None),
    #("C" , 1, None),
    #("Ca", 1, (3.5, 8)),
    #("Ca", 2, (3.5, 8)),
    #("Ce", 2, (-0.5, 4)),
    #("Co", 1, (0, 6)),
    #("Cr", 1, (0, 10)),
    #("Cu", 1, (1, 6)),
    #("Eu", 2, None),
    #("Fe", 1, (4, 8)),

    #("Fe", 2, (4, 8)),
    #("La", 2, (-2, +4)),
    #("Li", 1, None),
    #("Mg", 1, (5, 9)),
    #("Mn", 1, (0, 10)),
    #("Mn", 2, (0, 10)),
    #("Mo", 1, None),
    #("Na", 1, (2.5, 7.5)),
    #("Nd", 1, (-1, +5)),
    #("Nd", 2, (-1, +5)),
    #("Ni", 1, (3, 10)),
    #("Ni", 2, (3, 10)),
    #("O",  1, None),
    #("Pr", 2, (-3, 3)),
    #("Ru", 1, None),
    #("S",  1, None),
    #("Sc", 1, (0, 6)),
    #("Sc", 2, (0, 6)),
    #("Si", 1, None),
    #("Si", 2, None),
    #("Sm", 2, None),
    #("Sr", 1, None),
    #("Ti", 1, (0, 10)),
    ("Ti", 2, (0, 10)),
    ("V",  1, (1, 8)),
    ("V",  2, (1, 8)),
    ("Y",  2, (-2, +4)),
    ("Zn", 1, None),
    ("Zr", 1, (0, 6)),
    ("Zr", 2, (0, 6)),
    
]

"""
species = [
 #   ("Cr", 1, (0, 10)),
 #   ("Cr", 2, None),
 #   ("Ti", 1, None),
    ("Fe", 1, None),
]
"""
species = [
    ("Fe", 1, None)
]

species = [
    ("Si", 1, None),
    ("Si", 2, None),
]

for element, ion, absolute_extent in species:
    

    figures = OrderedDict([
        #("compare-bm", (ges.plot.benchmark_comparison,
        #    { "benchmark_filename": "benchmarks.yaml" })),
        #("compare-solar", ges.plot.solar_comparison),
        #("compare-m67-1194", ges.plot.m67_twin_comparison),
        #("benchmarks", (ges.plot.benchmark_line_abundances, { 
        #    "benchmark_filename": "benchmarks.yaml" })),
        #("differential-line-abundances", ges.plot.differential_line_abundances),
        #("differential-line-abundances-clipped", (
        #    ges.plot.differential_line_abundances, { "absolute_extent": absolute_extent })),
        #("abundance-heatmap", (ges.plot.transition_heatmap, {"column": "abundance"})),
        #("ew-heatmap", (ges.plot.transition_heatmap, {"column": "ew"})),
        #("mean-abundances", (ges.plot.mean_abundances, {
        #    "extent": absolute_extent })),

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
            "y_extent": absolute_extent,
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
            "y_extent": absolute_extent,
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
            "y_extent": absolute_extent,
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


        ("differential-line-abundances-wrt-snr", (
            ges.plot.differential_line_abundances_wrt_parameter, {
                "parameter": "snr",
                "x_extent": (0, 400) 
            }
        )),

        ("differential-line-abundances-wrt-snr-scaled", (
            ges.plot.differential_line_abundances_wrt_parameter, {
                "parameter": "snr",
                "x_extent": (0, 400),
                "scaled": True,
            }
        )),



        ("differential-line-abundances-wrt-teff", (
            ges.plot.differential_line_abundances_wrt_parameter, {
                "parameter": "teff",
                "x_extent": (3500, 7000)
            }
        )),
        ("differential-line-abundances-wrt-teff-scaled", (
            ges.plot.differential_line_abundances_wrt_parameter, {
                "parameter": "teff",
                "scaled": True,
                "x_extent": (3500, 7000)
            }
        )),

        ("differential-line-abundances-wrt-logg", (
            ges.plot.differential_line_abundances_wrt_parameter, {
                "parameter": "logg",
                "x_extent": (0, 5)
            }
        )),
        ("differential-line-abundances-wrt-logg-scaled", (
            ges.plot.differential_line_abundances_wrt_parameter, {
                "parameter": "logg",
                "scaled": True,
                "x_extent": (0, 5)
            }
        )),

        ("differential-line-abundances-wrt-feh", (
            ges.plot.differential_line_abundances_wrt_parameter, {
                "parameter": "feh",
                "x_extent": (-3, 0.5)
            }
        )),

        ("differential-line-abundances-wrt-feh-scaled", (
            ges.plot.differential_line_abundances_wrt_parameter, {
                "parameter": "feh",
                "scaled": True,
                "x_extent": (-3, 0.5)
            }
        )),

    ])


    """


    figures = OrderedDict([
           ("line-abundances-logx-wrt-feh-clip", (ges.plot.line_abundances, {
            "abundance_format": "log_x",
            "reference_column": "feh",
            "aux_column": "teff",
            "aux_extent": (3500, 7000),
            "y_extent": absolute_extent,
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
            "y_extent": absolute_extent,
            "highlight_flagged": True,
            "show_node_comparison": True,
            "show_line_comparison": True,
            "show_homogenised": True,
            })),




        ("homogenised-abundance-uncertainties", (
            ges.plot.homogenised_abundance_uncertainties, {
            })),
        ("line-abundances-logx-wrt-teff-clip", (ges.plot.line_abundances, {
            "abundance_format": "log_x",
            "reference_column": "teff",
            "aux_column": "logg",
            "aux_extent": (0, 5),
            "y_extent": absolute_extent,
            "highlight_flagged": True,
            "show_node_comparison": True,
            "show_line_comparison": True,
            "show_homogenised": True,
            })),
    ])


     ("compare-bm", (ges.plot.benchmark_comparison,
            { "benchmark_filename": "benchmarks.yaml" })),
        ("compare-solar", ges.plot.solar_comparison),
        ("compare-m67-1194", ges.plot.m67_twin_comparison),
        ("benchmarks", (ges.plot.benchmark_line_abundances, { 
            "benchmark_filename": "benchmarks.yaml" })),
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



        ("differential-line-abundances-wrt-snr", (
            ges.plot.differential_line_abundances_wrt_parameter, {
                "parameter": "snr",
                "x_extent": (0, 400) 
            }
        )),

        ("differential-line-abundances-wrt-snr-scaled", (
            ges.plot.differential_line_abundances_wrt_parameter, {
                "parameter": "snr",
                "x_extent": (0, 400),
                "scaled": True,
            }
        )),



        ("differential-line-abundances-wrt-teff", (
            ges.plot.differential_line_abundances_wrt_parameter, {
                "parameter": "teff",
                "x_extent": (3500, 7000)
            }
        )),
        ("differential-line-abundances-wrt-teff-scaled", (
            ges.plot.differential_line_abundances_wrt_parameter, {
                "parameter": "teff",
                "scaled": True,
                "x_extent": (3500, 7000)
            }
        )),

        ("differential-line-abundances-wrt-logg", (
            ges.plot.differential_line_abundances_wrt_parameter, {
                "parameter": "logg",
                "x_extent": (0, 5)
            }
        )),
        ("differential-line-abundances-wrt-logg-scaled", (
            ges.plot.differential_line_abundances_wrt_parameter, {
                "parameter": "logg",
                "scaled": True,
                "x_extent": (0, 5)
            }
        )),

        ("differential-line-abundances-wrt-feh", (
            ges.plot.differential_line_abundances_wrt_parameter, {
                "parameter": "feh",
                "x_extent": (-3, 0.5)
            }
        )),

        ("differential-line-abundances-wrt-feh-scaled", (
            ges.plot.differential_line_abundances_wrt_parameter, {
                "parameter": "feh",
                "scaled": True,
                "x_extent": (-3, 0.5)
            }
        )),
    """


    """

    figures = OrderedDict([
        ("line-abundances-logx-wrt-feh-clip", (ges.plot.line_abundances, {
            "abundance_format": "log_x",
            "reference_column": "feh",
            "aux_column": "teff",
            "aux_extent": (3500, 7000),
            "y_extent": absolute_extent,
            "highlight_flagged": True,
            "show_node_comparison": True,
            "show_line_comparison": True,
            "show_homogenised": True,
            })),

        ("homogenised-line-abundance-uncertainties", (
            ges.plot.homogenised_line_abundance_uncertainties, {
            })),
        ("homogenised-abundance-uncertainties", (
            ges.plot.homogenised_abundance_uncertainties, {
            })),

        ("compare-bm", (ges.plot.benchmark_comparison,
            { "benchmark_filename": "benchmarks.yaml" })),
        ("compare-solar", ges.plot.solar_comparison),
        ("compare-m67-1194", ges.plot.m67_twin_comparison),
        ("benchmarks", (ges.plot.benchmark_line_abundances, { 
            "benchmark_filename": "benchmarks.yaml" })),


    ])
    
    """

    """

        ("line-abundances-logx-wrt-logg-clip", (ges.plot.line_abundances, {
            "abundance_format": "log_x",
            "reference_column": "logg",
            "aux_column": "feh",
            "aux_extent": (-3, 1),
            "y_extent": absolute_extent,
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
            "y_extent": absolute_extent,
            "highlight_flagged": True,
            "show_node_comparison": True,
            "show_line_comparison": True,
            "show_homogenised": True,
            })),
    """

    figures = OrderedDict([
        ("differential-line-abundances", ges.plot.differential_line_abundances),
        ("differential-line-abundances-clipped", (
            ges.plot.differential_line_abundances, { "absolute_extent": absolute_extent })),
        ("scaled-differential-line-abundances", (
            ges.plot.differential_line_abundances, { "scaled": True })),
    ])

    try:
        make_figures(figures, element, ion)
    except:
        logger.exception("FAIL SNAIL")
        raise 

    plt.close("all")

