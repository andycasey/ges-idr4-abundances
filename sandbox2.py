
import os
import logging
from collections import OrderedDict
import matplotlib.pyplot as plt

import release

logger = logging.getLogger("ges")



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
                        figure.savefig(filename)
                        print("Created {}".format(filename))

            else:
                filename = os.path.join(directory, "{0}.{1}".format(
                    figure_name, format))
                fig.savefig(filename)
                print("Created {}".format(filename))


ges = release.DataRelease("arc")


species = [
    #("Si", 1, None),
    #("Si", 2, None),
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
        ("line-abundances-logx-wrt-teff", (ges.plot.line_abundances, {
            "reference_column": "teff",
            "abundance_format": "log_x",
            "aux_column": "logg",
            "aux_extent": (0, 5),
            "extent": None,
            "show_node_comparison": True,
            "show_line_comparison": True
        }))

    ])


    """
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
    """
    
    try:
        make_figures(figures, element, ion)
    except:
        logger.exception("FAIL SNAIL")
        raise 

    raise a
    
    plt.close("all")