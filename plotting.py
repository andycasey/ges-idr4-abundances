#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Plotting. """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import logging
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

import utils


#("compare-bm", (plot.compare_benchmarks,
#    { "benchmarks_filename": "benchmarks.yaml" })),
#("compare-solar", plot.compare_solar),
#("compare-m67-1194", plot.compare_m67_twin),

def latexify(label):
    """ A placeholder for a smart function to latexify common labels. """
    return label




def _compare_repeat_spectra(measurements, colors, homogenised_measurements=None,
    scaled=False, abundance_extent=None, bins=20, highlight_flagged=True,
    show_legend=True, x_column="SNR", reference_abundance=None, 
    reference_uncertainty=None, reference_label=None):

    nodes = sorted(set(measurements["node"]))
    wavelengths = sorted(set(measurements["wavelength"]))
    N_nodes, N_wavelengths = len(nodes), len(wavelengths)

    column = "scaled_abundance" if scaled else "abundance"

    # We should be getting the colours passed to us.
    bin_min, bin_max = abundance_extent \
        or (np.nanmin(measurements[column]), np.nanmax(measurements[column]))

    if reference_abundance is not None:
        bin_min = min([bin_min, reference_abundance - 0.1])
        bin_max = max([bin_max, reference_abundance + 0.1])

    if bin_min == bin_max:
        bin_min, bin_max = (bin_min - 0.5, bin_max + 0.5)

    if homogenised_measurements is not None:
        bin_min = min([bin_min, np.nanmin(homogenised_measurements["abundance"])])
        bin_max = max([bin_max, np.nanmax(homogenised_measurements["abundance"])])

    hist_kwds = {
        "histtype": "step",
        "bins": np.linspace(bin_min, bin_max, bins + 1),
        "normed": True,
        "lw": 2
    }
    scatter_kwds = {
        "s": 50,
        "zorder": 10,
    }
    span_kwds = {
        "alpha": 0.5,
        "edgecolor": "none",
        "zorder": -100,
    }
    hline_kwds = { "zorder": -10, "lw": 2}

    # Histograms be square as shit.
    K = 2
    scale = 2.0
    wspace, hspace = 0.45, 0.3
    lb, tr = 0.5, 0.2
    ys = scale * N_wavelengths + scale * (N_wavelengths - 1) * wspace
    xs = scale * K + scale * (K - 1) * hspace
    
    xdim = lb * scale + xs + tr * scale
    ydim = lb * scale + ys + tr * scale

    fig, axes = plt.subplots(N_wavelengths, K, figsize=(xdim, ydim))
    fig.subplots_adjust(
        left=(lb * scale)/xdim,
        bottom=(lb * scale)/ydim,
        right=(lb * scale + xs)/xdim,
        top=(tr * scale + ys)/ydim,
        wspace=wspace, hspace=hspace)

    flagged = measurements["flags"] > 0
    for i, (wavelength, (ax_hist, ax_snr)) in enumerate(zip(wavelengths, axes)):

        # Show the normalised distribution of absolute abundances from all nodes
        ok = (measurements["wavelength"] == wavelength) * \
            np.isfinite(measurements[column])
        ax_hist.set_title(latexify(wavelength))

        # For each node, show the histogram of absolute abundances.
        for j, node in enumerate(nodes):
            node_match = ok * (measurements["node"] == node)
            if np.any(np.isfinite(measurements[column][node_match])):
                ax_hist.hist(measurements[column][node_match],
                    color=colors[node], **hist_kwds)

            # For the legend.
            ax_hist.plot([], [], color=colors[node], label=node)

        # For each node, show the points and uncertainties as a function of S/N
        for j, node in enumerate(nodes):
            node_match = ok * (measurements["node"] == node)

            # Any flagged lines?
            if highlight_flagged and any(flagged * node_match):
                # Show them differently.
                match = flagged * node_match
                x = measurements["snr"][match]
                y = measurements[column][match]
                ax_snr.errorbar(x, y, yerr=measurements["e_abundance"][match],
                    lc="k", ecolor="k", aa=True, fmt=None, mec="k", mfc="w",
                    ms=6, zorder=1)
                kwds = scatter_kwds.copy()
                kwds["edgecolor"] = "r"
                kwds["lw"] = 2
                ax_snr.scatter(x, y, facecolor=colors[node], **kwds)

                node_match *= ~flagged
    
            x = measurements["snr"][node_match]
            y = measurements[column][node_match]
            ax_snr.errorbar(x, y, yerr=measurements["e_abundance"][node_match],
                lc="k", ecolor="k", aa=True, fmt=None, mec="k", mfc="w", ms=6,
                zorder=1)
            ax_snr.scatter(x, y, facecolor=colors[node], **scatter_kwds)

            # Show some mean + std.dev range.
            mean, sigma = np.nanmean(y), np.nanstd(y)
            ax_snr.axhline(mean, c=colors[node], **hline_kwds)
            if np.isfinite(mean * sigma):
                ax_snr.axhspan(mean - sigma, mean + sigma,
                    facecolor=colors[node], **span_kwds)

            ax_snr.set_ylim(ax_hist.get_xlim())

        if homogenised_measurements is not None:
            homogenised_match = \
                (homogenised_measurements["wavelength"] == wavelength) * \
                np.isfinite(homogenised_measurements["abundance"])

            x = homogenised_measurements["snr"][homogenised_match]
            y = homogenised_measurements["abundance"][homogenised_match]
            yerr = homogenised_measurements["e_abundance"][homogenised_match]

            kwds = scatter_kwds.copy()
            kwds.update({
                "facecolor": "w",
                "zorder": 100,
                "s": 50,
                "lw": 2
            })
            ax_snr.scatter(x, y, label="Homogenised", **kwds)
            ax_snr.errorbar(x, y, yerr=yerr, lc="k", ecolor="k", aa=True,
                fmt=None, mec="k", mfc="w", ms=6, lw=2, zorder=99)

        ax_hist.xaxis.set_major_locator(MaxNLocator(5))
        ax_hist.yaxis.set_major_locator(MaxNLocator(5))
        ax_hist.set_ylim(0, ax_hist.get_ylim()[1])

        label = r"{0} {1}".format(measurements["element"][0],
            measurements["ion"][0])
        ax_snr.set_ylabel(latexify(label))
        ax_hist.set_ylabel(latexify(label))
        
        if not ax_hist.is_last_row():
            ax_hist.set_xticklabels([])
            ax_snr.set_xticklabels([])

        else:
            ax_snr.set_xlabel(latexify(x_column))
            ax_hist.set_xlabel(latexify(label))

        if ax_hist.is_first_row():
            ax_hist.legend(loc="upper left", frameon=False, fontsize=12)

        ax_snr.xaxis.set_major_locator(MaxNLocator(5))
        ax_snr.yaxis.set_major_locator(MaxNLocator(5))

        # Show the reference abundance.
        if reference_abundance is not None:
            ax_snr.axhline(reference_abundance, c="k", lw=3,
                label=reference_label, zorder=1)
            if reference_uncertainty is not None:
                ax_snr.axhspan(reference_abundance - reference_uncertainty,
                    reference_abundance + reference_uncertainty,
                    facecolor="#666666", alpha=0.5, zorder=-100,
                    edgecolor="#666666")
            if ax_snr.is_first_row() and reference_label is not None \
            and show_legend:
                ax_snr.legend(loc="upper right", frameon=False, fontsize=12)

    # Set S/N axes on the same x-scale.
    xlims = np.array([ax.get_xlim() for ax in axes[:, 1]])
    [ax.set_xlim(xlims[:, 0].min(), xlims[:,1].max()) for ax in axes[:, 1]]

    return fig



class AbundancePlotting(object):

    def __init__(self, release):
        self.release = release

    def solar_comparison(self, element, ion, scaled=False, 
        highlight_flagged=True, show_homogenised=True, bins=20,
        abundance_extent=None, **kwargs):
        """
        Show distributions of abundances for each node and each line for the Sun

        :param element:
            The element name to show comparisons for.

        :type element:
            str

        :param ion:
            The ionisation stage of the element to show (1 = neutral).

        :type ion:
            int

        :param scaled: [optional]
            Show scaled abundances. Alternative is to show unscaled abundances.

        :type scaled:
            bool

        :param highlight_flagged: [optional]
            Show a red outline around measurements that are flagged.

        :type highlight_flagged:
            bool

        :param show_homogenised: [optional]
            Show homogenised line abundances.

        :type show_homogenised:
            bool

        :param bins: [optional]
            The number of bins to show in the histograms.

        :type bins:
            int

        :param abundance_extent: [optional]
            The lower and upper range to show in abundances.

        :type abundance_extent:
            two length tuple of floats
        """

        measurements = self.release.retrieve_table(
            """SELECT * FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) 
            cname, snr FROM node_results ORDER BY cname) n ON (l.element = '{0}'
            AND l.ion = '{1}' AND l.cname = n.cname AND l.cname LIKE 'ssssss%')
            """.format(element, ion))
        if measurements is None: return None

        homogenised_measurements = self.release.retrieve_table(
            """SELECT * FROM homogenised_line_abundances l JOIN (SELECT
            DISTINCT ON (cname) cname, snr FROM node_results ORDER BY cname)
            n ON (l.element = '{0}' AND l.ion = '{1}' AND l.cname = n.cname
            AND l.cname LIKE 'ssssss%')""".format(element, ion)) \
            if show_homogenised else None

        return _compare_repeat_spectra(measurements, self.release.node_colors,
            homogenised_measurements, scaled=scaled, x_column="SNR",
            abundance_extent=abundance_extent, bins=bins,
            highlight_flagged=highlight_flagged,
            reference_abundance=utils.solar_abundance(element),
            reference_uncertainty=None, reference_label="Solar", **kwargs)


    def m67_twin_comparison(self, element, ion, scaled=False, 
        highlight_flagged=True, show_homogenised=True, bins=20,
        abundance_extent=None, **kwargs):
        """
        Show distributions of abundances for each node and each line for the 
        solar twin in M67 (M67-1194).

        :param element:
            The element name to show comparisons for.

        :type element:
            str

        :param ion:
            The ionisation stage of the element to show (1 = neutral).

        :type ion:
            int

        :param scaled: [optional]
            Show scaled abundances. Alternative is to show unscaled abundances.

        :type scaled:
            bool

        :param highlight_flagged: [optional]
            Show a red outline around measurements that are flagged.

        :type highlight_flagged:
            bool

        :param show_homogenised: [optional]
            Show homogenised line abundances.

        :type show_homogenised:
            bool

        :param bins: [optional]
            The number of bins to show in the histograms.

        :type bins:
            int

        :param abundance_extent: [optional]
            The lower and upper range to show in abundances.

        :type abundance_extent:
            two length tuple of floats
        """

        measurements = self.release.retrieve_table(
            """SELECT * FROM line_abundances l JOIN (SELECT DISTINCT ON (cname) 
            cname, snr, ges_fld, object FROM node_results ORDER BY cname) n ON
            (l.element = '{0}' AND l.ion = '{1}' AND l.cname = n.cname 
                AND n.ges_fld LIKE 'M67%' AND n.object = '1194')"""\
            .format(element, ion))
        if measurements is None: return None

        homogenised_measurements = self.release.retrieve_table(
            """SELECT * FROM homogenised_line_abundances l JOIN (SELECT
            DISTINCT ON (cname) cname, snr, ges_fld, object FROM node_results
            ORDER BY cname) n ON (l.element = '{0}' AND l.ion = '{1}' AND 
            l.cname = n.cname AND n.ges_fld LIKE 'M67%' AND n.object = '1194')
            """.format(element, ion)) if show_homogenised else None

        return _compare_repeat_spectra(measurements, self.release.node_colors,
            homogenised_measurements, scaled=scaled, x_column="SNR",
            abundance_extent=abundance_extent, bins=bins,
            highlight_flagged=highlight_flagged,
            reference_abundance=utils.solar_abundance(element),
            reference_uncertainty=None, reference_label="Solar", **kwargs)


    def differential_line_abundances(self, element, ion, scaled=False, bins=50,
        absolute_abundance_extent=None, differential_abundance_extent=(-0.5, 0.5),
        ignore_flagged=True, show_legend=True, **kwargs):
        """
        Show histograms of the absolute and differential line abundances for a
        given element and ion.
        
        :param element:
            The atomic element of interest.

        :type element:
            str

        :param ion:
            The ionisation stage of the species of interest (1 indicates neutral)

        :type ion:
            int

        :param scaled: [optional]
            Show scaled abundances. Alternative is to show unscaled abundances.

        :type scaled:
            bool

        :param bins: [optional]
            The number of bins to show in the histograms.

        :type bins:
            int

        :param absolute_abundance_extent: [optional]
            The lower and upper range to show in absolute abundances.

        :type absolute_abundance_extent:
            None or two length tuple of floats

        :param differential_abundance_extent: [optional]
            The lower and upper range to show in differential abundances.

        :type differential_abundance_extent:
            None or a two-length tuple of floats

        :param ignore_flagged: [optional]
            Do not include measurements that have been flagged as poor quality.
        
        :type ignore_flagged:
            bool
        """

        data = self.release.retrieve_table("""SELECT * FROM line_abundances
            WHERE element = %s AND ion = %s {flag_query} ORDER BY wavelength ASC
            """.format(flag_query="AND flags = 0" if ignore_flagged else ""),
            (element, ion))
        if data is None: return

        # Use only finite measurements.
        column = "scaled_abundance" if scaled else "abundance"
        use = np.isfinite(data[column]) * (data["upper_abundance"] == 0)
        if not any(use): return None

        data = data[use]
        wavelengths = set(data["wavelength"])

        K, N_lines = 3, len(wavelengths)

        # Histograms be square as shit.
        scale = 2.0
        wspace, hspace = 0.3, 0.2
        lb, tr = 0.5, 0.2
        ys = scale * N_lines + scale * (N_lines - 1) * wspace
        xs = scale * K + scale * (K - 1) * hspace
        
        xdim = lb * scale + xs + tr * scale
        ydim = lb * scale + ys + tr * scale

        fig, axes = plt.subplots(N_lines, K, figsize=(xdim, ydim))
        fig.subplots_adjust(
            left=(lb * scale)/xdim,
            bottom=(lb * scale)/ydim,
            right=(lb * scale + xs)/xdim,
            top=(tr * scale + ys)/ydim,
            wspace=wspace, hspace=hspace)


        bin_min, bin_max = absolute_abundance_extent \
            or (data[column].min(), data[column].max())
       
        if abs(bin_min - bin_max) < 0.005:
            bin_min, bin_max = bin_min - 0.5, bin_min + 0.5

        hist_kwds = {
            "histtype": "step",
            "bins": np.linspace(bin_min, bin_max, bins + 1),
            "normed": True,
        }
        full_distribution_color = "k"
        comp_distribution_color = "#666666"

        for i, (ax, wavelength) in enumerate(zip(axes.T[0], wavelengths)):
            # Show the full distribution.
            ax.hist(data[column], color=full_distribution_color, **hist_kwds)

            # Show the distribution for this wavelength.
            match = data["wavelength"] == wavelength
            ax.hist(data[column][match], color=comp_distribution_color,
                **hist_kwds)

            ax.xaxis.set_major_locator(MaxNLocator(5))
            ax.yaxis.set_major_locator(MaxNLocator(5))
            if not ax.is_last_row():
                ax.set_xticklabels([])

            else:
                ax.set_xlabel(latexify("{0} {1}".format(element, ion)))
            
            ax.text(0.95, 0.95, latexify(len(data)), transform=ax.transAxes,
                verticalalignment="top", horizontalalignment="right",
                color=full_distribution_color)
            ax.text(0.95, 0.95, latexify("\n{}".format(match.sum())),
                transform=ax.transAxes, verticalalignment="top",
                horizontalalignment="right",
                color=comp_distribution_color)

        # Plot the full distribution of differential abundances.
        # (but first,..) determine the full matrix of differential abundances.
        X, nodes, diff_data = self.release._match_species_abundances(
            element, ion, scaled=scaled, include_flagged_lines=not ignore_flagged)

        # Calculate the full differential abundances.
        X_diff, indices = utils.calculate_differential_abundances(X,
            full_output=True)
        X_diff = X_diff[np.isfinite(X_diff)]

        b_min, b_max = differential_abundance_extent \
            or (np.nanmin(X_diff), np.nanmax(X_diff))
        hist_kwds["bins"] = np.linspace(b_min, b_max, bins + 1)

        for i, (ax, wavelength) in enumerate(zip(axes.T[1], wavelengths)):
            if ax.is_first_row():
                ax.text(0.05, 0.95,
                    r"$\mu = {0:.2f}$" "\n" r"$\sigma = {1:.2f}$".format(
                    np.nanmean(X_diff), np.nanstd(X_diff)),  fontsize=10,
                    transform=ax.transAxes, color=full_distribution_color,
                    verticalalignment="top", horizontalalignment="left")
                    
            # Show the full distribution of differential abundances.
            ax.hist(X_diff, color=full_distribution_color, **hist_kwds)

            # Show the distribution of differential abundances for this wavelength.
            match = diff_data["wavelength"] == wavelength
            X_diff_wavelength = utils.calculate_differential_abundances(X[match],
                full_output=False).flatten()
            if np.isfinite(X_diff_wavelength).sum() > 0:
                ax.hist(X_diff_wavelength, color=comp_distribution_color,
                    **hist_kwds)

            ax.set_title(latexify(wavelength))
            ax.xaxis.set_major_locator(MaxNLocator(5))
            ax.yaxis.set_major_locator(MaxNLocator(5))
            if not ax.is_last_row():
                ax.set_xticklabels([])
            else:
                ax.set_xlabel(r"$\Delta${0} {1}".format(element, ion))

            ax.text(0.95, 0.95, latexify(X_diff.size), transform=ax.transAxes,
                verticalalignment="top", horizontalalignment="right",
                color=full_distribution_color)
            ax.text(0.95, 0.95, 
                latexify("\n{}".format(np.isfinite(X_diff_wavelength).sum())),
                verticalalignment="top", horizontalalignment="right",
                color=comp_distribution_color, transform=ax.transAxes)
            
        # Plot the node-to-node distribution of the differential abundances.
        # Node X compared to all others
        # Node Y compared to all others, etc.
        colors = self.release.node_colors
        for i, (ax, wavelength) in enumerate(zip(axes.T[2], wavelengths)):
            # Show the full distribution of differential abundances.
            # Show the distribution of differential abundances for each node.
            match = diff_data["wavelength"] == wavelength
            X_diff_wavelength = utils.calculate_differential_abundances(X[match],
                full_output=False)

            if np.any(np.isfinite(X_diff_wavelength)):
                ax.hist(X_diff_wavelength.flatten(),
                    color=full_distribution_color, **hist_kwds)

            else:
                ax.set_ylim(0, 1) # For prettyness

            ax.text(0.95, 0.95, np.isfinite(X_diff_wavelength).sum(),
                transform=ax.transAxes, color=full_distribution_color,
                verticalalignment="top", horizontalalignment="right")
            

            for j, node in enumerate(nodes):
                ax.plot([], [], label=node, c=colors[node])

                # This -1, +1 business ensures the third column always contains
                # Node - someone else.
                X_diff_wavelength_node = np.hstack([[-1, +1][j == idx[0]] * \
                    X_diff_wavelength[:, k].flatten() \
                    for k, idx in enumerate(indices) if j in idx])

                if np.any(np.isfinite(X_diff_wavelength_node)):
                    ax.hist(X_diff_wavelength_node, color=colors[node],
                        **hist_kwds)
                
                ax.text(0.95, 0.95, latexify("{0}{1}".format("\n"*(j + 1), 
                        np.isfinite(X_diff_wavelength_node).sum())),
                    transform=ax.transAxes, verticalalignment="top",
                    horizontalalignment="right", color=colors[node])

            if np.any(np.isfinite(X_diff_wavelength)):
                ax.text(0.05, 0.95, 
                    r"$\mu = {0:.2f}$" "\n" r"$\sigma = {1:.2f}$".format(
                    np.nanmean(X_diff_wavelength), np.nanstd(X_diff_wavelength)),
                    fontsize=10, transform=ax.transAxes, 
                    color=full_distribution_color,
                    verticalalignment="top", horizontalalignment="left")

            ax.xaxis.set_major_locator(MaxNLocator(5))
            ax.yaxis.set_major_locator(MaxNLocator(5))
            if not ax.is_last_row():
                ax.set_xticklabels([])
            else:
                ax.set_xlabel(r"$\Delta${0} {1}".format(element, ion))

        if show_legend:
            axes.T[0][0].legend(*axes.T[2][0].get_legend_handles_labels(),
                loc="upper left", frameon=False, fontsize=10)

        return fig