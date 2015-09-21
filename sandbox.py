#!/usr/bin/python

""" Make plots for an element. """

import os
import psycopg2 as pg
import plot

import matplotlib.pyplot as plt



def make_figures(figures, database, element, ion, format="png"):

    directory = "figures/{0}{1}".format(element.upper(), ion)
    if not os.path.exists(directory): os.mkdir(directory)
    for figure_name, details in figures.items():
        
        if not isinstance(details, (tuple, )):
            details = (details, )

        command = details[0]
        kwargs = {} if len(details) == 1 else details[1]
        

        print("Doing {0} with {1}".format(figure_name, kwargs))
        fig = command(database, element, ion, **kwargs)
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


if __name__ == "__main__":

    absolute_extent = (5, 9)
    figures = {
        "compare-bm": (plot.compare_benchmarks,
            { "benchmarks_filename": "benchmarks.yaml" }),
        "compare-solar": plot.compare_solar,
        "compare-m67-1194": plot.compare_m67_twin,
        "percentile": plot.percentiles,
        "differential-line-abundances": plot.differential_line_abundances,
        "differential-line-abundances-clipped": (
            plot.differential_line_abundances, { "absolute_extent": absolute_extent }),     
        "line-abundances-logx-wrt-teff": (plot.line_abundances, {
            "reference_column": "teff",
            "aux_column": "logg",
            "aux_extent": (3500, 7500),
            "extent": None,
            "show_node_comparison": False,
            "show_line_comparison": False,
            }),
        "line-abundances-logx-wrt-logg": (plot.line_abundances, {
            "reference_column": "logg",
            "aux_column": "feh",
            "aux_extent": (0, 5),
            "extent": None,
            "show_node_comparison": False,
            "show_line_comparison": False,
            }),
        "line-abundances-logx-wrt-feh": (plot.line_abundances, {
            "reference_column": "feh",
            "aux_column": "teff",
            "aux_extent": (-3, 1),
            "extent": None,
            "show_node_comparison": False,
            "show_line_comparison": False,
            }),

        "line-abundances-xfe-wrt-teff": (plot.line_abundances, {
            "reference_column": "teff", "aux_column": "logg",
            "aux_extent": (0, 5) }),
        "line-abundances-xfe-wrt-logg": (plot.line_abundances, {
            "reference_column": "logg", "aux_column": "feh",
            "aux_extent": (-3, 1) }),
        "line-abundances-xfe-wrt-feh": (plot.line_abundances, {
            "reference_column": "feh", "aux_column": "teff",
            "aux_extent": (3500, 7000) }),


        "abundance-heatmap": (plot.transition_heatmap, {"column": "abundance"}),
        "ew-heatmap": (plot.transition_heatmap, {"column": "ew"}),
        #"abundance-covariance": (plot.transition_covariance, {"column": "abundance"}),
        #"ew-covariance": (plot.transition_covariance, {"column": "ew"}),
        "mean-abundance-sp": plot.mean_abundance_against_stellar_parameters,
        "mean-abundance-differences": (plot.mean_abundance_differences, {
            "extent": absolute_extent }),
        "line-abundances-rew": (plot.all_node_individual_line_abundance_differences,
            {"rew_on_x_axis": True}),
        "line-abundances": (plot.all_node_individual_line_abundance_differences,
            {"rew_on_x_axis": False}),
        "line-abundances-rew-clip": (plot.all_node_individual_line_abundance_differences,
            {"rew_on_x_axis": True, "x_extent": (-6.5, -4), "y_extent": (-1.5, 1.5), "vmin": 6, "vmax": 9}),
        "line-abundances-clip": (plot.all_node_individual_line_abundance_differences,
            {"rew_on_x_axis": False, "x_extent": absolute_extent, "y_extent": (-1.5, 1.5)}),    
    }

    species = [
        ("Si", 1), # Done
        ("Si", 2), # Done
    ]
    """
        ("Al", 1),
        ("Al", 3),
        ("Na", 1),
        ("S", 1),
        ("S", 2),
        ("S", 3),
        ("Ne", 1),
        ("Ne", 2),
        
        ("Mg", 1),
        ("Ca", 1),
        ("Ca", 2),

        ("Ni", 1),
        ("Zn", 1),
        ("Cr", 1),
        ("Cr", 2),
        ("Sc", 1),
        ("Sc", 2),
        ("Cu", 1),
        ("Co", 1),
        ("Mn", 1),
        ("Zr", 1),
        ("Zr", 2),
        ("V", 1),
        ("V", 2),
        ("O", 1),
        ("O", 2),
        ("Ti", 1),
        ("Ti", 2),
        ("Fe", 1),
        ("Fe", 2),
        ("Fe", 3),
        ("La", 2),
        ("Ba", 2),
        ("Y", 1),
        ("Y", 2),
        ("C", 1),
        ("C", 2),
        ("C", 3),
        #("C_C", )
        #("N_C", )
        ("N", 1),
        ("N", 2),
        ("Li", 1),
        ("Ru", 1),
        ("Sr", 1),
        ("Dy", 2),
        ("Ce", 2),
        ("Pr", 2),
        ("Nb", 1),
        ("Nd", 2),
        ("Mo", 1),
        ("Eu", 2),
        ("Gd", 2),
        ("Sm", 2),
    ]
    """

    database = pg.connect(dbname="arc")


    #make_figures(figures, database, "Si", 1)

    for element, ion in species:
        try:
            make_figures(figures, database, element, ion)
        except:
            raise
            None
        plt.close("all")

    #for i in range(1, 5):
    #    make_figures(figures, database, "Si", i)