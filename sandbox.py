#!/usr/bin/python

""" Make plots for an element. """

import os
import psycopg2 as pg
import plot




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


    figures = {
        "abundance-heatmap": (plot.transition_heatmap, {"column": "abundance"}),
        "ew-heatmap": (plot.transition_heatmap, {"column": "ew"}),
        "abundance-covariance": (plot.transition_covariance, {"column": "abundance"}),
        "ew-covariance": (plot.transition_covariance, {"column": "ew"}),
        "mean-abundance-sp": plot.mean_abundance_against_stellar_parameters,
        "mean-abundance-differences": plot.mean_abundance_differences,
        "line-abundances-rew": (plot.all_node_individual_line_abundance_differences,
            {"rew_on_x_axis": True}),
        "line-abundances": (plot.all_node_individual_line_abundance_differences,
            {"rew_on_x_axis": False}),
        "line-abundances-rew-clip": (plot.all_node_individual_line_abundance_differences,
            {"rew_on_x_axis": True, "x_extent": (-6.5, -4), "y_extent": (-1.5, 1.5), "vmin": 6, "vmax": 9}),
        "line-abundances-clip": (plot.all_node_individual_line_abundance_differences,
            {"rew_on_x_axis": False, "x_extent": (6, 9), "y_extent": (-1.5, 1.5)}),    
    }

    database = pg.connect(dbname="arc")
    for i in range(3, 5):
        make_figures(figures, database, "Si", i)