#!/usr/bin/python

""" Make plots for an element. """

import os
import psycopg2 as pg
import plot

OUTPUT = "figures/"

def element_figures(nodes, element):

    directory = os.path.join(OUTPUT, element)
    if not os.path.exists(directory):
        os.mkdir(directory)
    
    # Make node-specific figures for this element.
    for node in nodes:
        fig = plot.node_element_abundance_against_stellar_parameters(
            node, element)
        filename = os.path.join(directory, "{0}-stellar-parameters.png".format(
            node[0].header["NODE1"]))
        fig.savefig(filename)
        print("Created {}".format(filename))


if __name__ == "__main__":

    format = "png"
    element, ion = "Si", 1
    directory = "figures/{0}{1}".format(element.upper(), ion)


    figures = {
        "abundance-heatmap": (plot.transition_heatmap, {"column": "abundance"}),
        "ew-heatmap": (plot.transition_heatmap, {"column": "ew"}),
        "abundance-covariance": (plot.transition_covariance, {"column": "abundance"}),
        "ew-covariance": (plot.transition_covariance, {"column": "ew"}),
        "mean-abundance-sp": plot.mean_abundance_against_stellar_parameters,
    }
    database = pg.connect(dbname="arc")
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
                    figure.savefig(filename)
                    print("Created {}".format(filename))

            else:
                filename = os.path.join(directory, "{0}.{1}".format(
                    figure_name, format))
                fig.savefig(filename)
                print("Created {}".format(filename))
