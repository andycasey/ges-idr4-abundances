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

    REMAKE = False

    from collections import OrderedDict
    

    species = [
        ("Si", 1, (5, 10)), # Done
        ("Si", 2, (5, 10)), # Done
        ("Al", 1, (3, 7.5)),
        ("Al", 3, (3, 7.5)),
        ("Na", 1, None),
        ("S", 1, None),
        ("S", 2, None),
        ("S", 3, None),
        ("Ne", 1, None),
        ("Ne", 2, None),
        
        ("Mg", 1, None),
        ("Ca", 1, None),
        ("Ca", 2, None),

        ("Ni", 1, None),
        ("Zn", 1, None),
        ("Cr", 1, None),
        ("Cr", 2, None),
        ("Sc", 1, None),
        ("Sc", 2, None),
        ("Cu", 1, None),
        ("Co", 1, None),
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
        ("Sm", 2, None),
    ]

    database = pg.connect(dbname="arc")


    #make_figures(figures, database, "Si", 1)

    for element, ion, absolute_extent in species:

        figures = OrderedDict([
            ("compare-bm", (plot.compare_benchmarks,
                { "benchmarks_filename": "benchmarks.yaml" })),
            ("compare-solar", plot.compare_solar),
            ("compare-m67-1194", plot.compare_m67_twin),
            ("benchmarks", (plot.benchmark_line_abundances, {
                "benchmark_filename": "benchmarks.yaml" })),
            ("percentile", plot.percentiles),
            ("differential-line-abundances", plot.differential_line_abundances),
            ("differential-line-abundances-clipped", (
                plot.differential_line_abundances, { "absolute_extent": absolute_extent })),

            ("differential-line-abundances-wrt-teff", (
                plot.differential_line_abundances_wrt_x,
                { "parameter": "teff", "x_extent": (3500, 7000) })),
            ("differential-line-abundances-wrt-logg", (
                plot.differential_line_abundances_wrt_x,
                { "parameter": "logg", "x_extent": (0, 5) })),
            ("differential-line-abundances-wrt-feh", (
                plot.differential_line_abundances_wrt_x,
                { "parameter": "feh", "x_extent": (-3, 0.5) })),
            ("line-abundances-logx-wrt-teff", (plot.line_abundances, {
                "reference_column": "teff",
                "abundance_format": "log_x",
                "aux_column": "logg",
                "aux_extent": (0, 5),
                "extent": None,
                "show_node_comparison": False,
                "show_line_comparison": False,
                })),
            ("line-abundances-logx-wrt-logg", (plot.line_abundances, {
                "abundance_format": "log_x",
                "reference_column": "logg",
                "aux_column": "feh",
                "aux_extent": (-3, 1),
                "extent": None,
                "show_node_comparison": False,
                "show_line_comparison": False,
                })),
            ("line-abundances-logx-wrt-feh", (plot.line_abundances, {
                "abundance_format": "log_x",
                "reference_column": "feh",
                "aux_column": "teff",
                "aux_extent": (3500, 7000),
                "extent": None,
                "show_node_comparison": False,
                "show_line_comparison": False,
                })),
            ("line-abundances-logx-wrt-teff-clip", (plot.line_abundances, {
                "abundance_format": "log_x",
                "reference_column": "teff",
                "aux_column": "logg",
                "aux_extent": (0, 5),
                "extent": absolute_extent,
                "show_node_comparison": False,
                "show_line_comparison": False,
                })),
            ("line-abundances-logx-wrt-logg-clip", (plot.line_abundances, {
                "abundance_format": "log_x",
                "reference_column": "logg",
                "aux_column": "feh",
                "aux_extent": (-3, 1),
                "extent": absolute_extent,
                "show_node_comparison": False,
                "show_line_comparison": False,
                })),
            ("line-abundances-logx-wrt-feh-clip", (plot.line_abundances, {
                "abundance_format": "log_x",
                "reference_column": "feh",
                "aux_column": "teff",
                "aux_extent": (3500, 7000),
                "extent": absolute_extent,
                "show_node_comparison": False,
                "show_line_comparison": False,
                })),
            ("line-abundances-xfe-wrt-teff", (plot.line_abundances, {
                "abundance_format": "x_fe",
                "reference_column": "teff", "aux_column": "logg",
                "aux_extent": (0, 5) })),
            ("line-abundances-xfe-wrt-logg", (plot.line_abundances, {
                "abundance_format": "x_fe",
                "reference_column": "logg", "aux_column": "feh",
                "aux_extent": (-3, 1) })),
            ("line-abundances-xfe-wrt-feh", (plot.line_abundances, {
                "abundance_format": "x_fe",
                "reference_column": "feh", "aux_column": "teff",
                "aux_extent": (3500, 7000) })),

            ("abundance-heatmap", (plot.transition_heatmap, {"column": "abundance"})),
            ("ew-heatmap", (plot.transition_heatmap, {"column": "ew"})),
            #"abundance-covariance": (plot.transition_covariance, {"column": "abundance"}),
            #"ew-covariance": (plot.transition_covariance, {"column": "ew"}),
            ("mean-abundance-sp", plot.mean_abundance_against_stellar_parameters),
            ("mean-abundance-differences", (plot.mean_abundance_differences, {
                "extent": absolute_extent })),
            ("line-abundances-rew", (plot.all_node_individual_line_abundance_differences,
                {"rew_on_x_axis": True})),
            ("line-abundances", (plot.all_node_individual_line_abundance_differences,
                {"rew_on_x_axis": False})),
            ("line-abundances-rew-clip", (plot.all_node_individual_line_abundance_differences,
                {"rew_on_x_axis": True, "x_extent": (-6.5, -4), "y_extent": (-1.5, 1.5), "vmin": 6, "vmax": 9})),
            ("line-abundances-clip", (plot.all_node_individual_line_abundance_differences,
                {"rew_on_x_axis": False, "x_extent": absolute_extent, "y_extent": (-1.5, 1.5)})),    
        ])

        
        figures = OrderedDict([
            ("differential-line-abundances-wrt-teff", (
                plot.differential_line_abundances_wrt_x,
                { "parameter": "teff", "x_extent": (3500, 7000) })),
            ("differential-line-abundances-wrt-logg", (
                plot.differential_line_abundances_wrt_x,
                { "parameter": "logg", "x_extent": (0, 5) })),
            ("differential-line-abundances-wrt-feh", (
                plot.differential_line_abundances_wrt_x,
                { "parameter": "feh", "x_extent": (-3, 0.5) })),

        ])


        figures = OrderedDict([
            ("differential-line-abundances", plot.differential_line_abundances),
        ])

        try:
            make_figures(figures, database, element, ion)
        except:
            raise
            None
        plt.close("all")

        # TODO: Check that we didn't miss anything?

    #for i in range(1, 5):
    #    make_figures(figures, database, "Si", i)