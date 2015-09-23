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
