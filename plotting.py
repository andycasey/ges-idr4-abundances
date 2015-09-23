#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Plotting. """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import logging
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from matplotlib.ticker import MaxNLocator

import colormaps
import utils

logger = logging.getLogger("ges")


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

        :param show_legend: [optional]
            Show a legend with the node colours.

        :type show_legend:
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
        X, nodes, diff_data = self.release._match_species_abundances(element,
            ion, scaled=scaled, include_flagged_lines=not ignore_flagged)

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


    def benchmark_comparison(self, element, ion, benchmark_filename, 
        scaled=False, highlight_flagged=True, show_homogenised=True, bins=20,
        abundance_extent=None, **kwargs):
        """
        Show distributions of abundances for each node and each line for all of
        the benchmarks (in separate figures).

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

        :returns:
            Returns a dictionary containing the figures as values.
        """

        figures = {}
        for benchmark in utils.load_benchmarks(benchmark_filename):

            if element not in benchmark.abundances:
                logger.warn("No elemental abundance for {0} in {1}".format(
                    element, benchmark.name))
                figures[benchmark.name] = None
                continue

            measurements = self.release.retrieve_table(
                """SELECT * FROM line_abundances l JOIN (SELECT DISTINCT ON
                (cname) cname, snr, ges_fld, ges_type, object FROM node_results
                ORDER BY cname) n ON (l.element = '{0}' AND l.ion = '{1}' AND
                l.cname = n.cname AND ges_type LIKE '%_BM' AND
                (ges_fld ILIKE '{2}%' OR n.object ILIKE '{2}%'))""".format(
                    element, ion, benchmark.name))
            
            homogenised_measurements = self.release.retrieve_table(
                """SELECT * FROM homogenised_line_abundances l 
                JOIN (SELECT DISTINCT ON (cname) cname, snr, ges_fld, ges_type,
                object FROM node_results ORDER BY cname) n
                ON (l.element = '{0}' AND l.ion = '{1}' AND
                l.cname = n.cname AND ges_type LIKE '%_BM' AND
                (ges_fld ILIKE '{2}%' OR n.object ILIKE '{2}%'))""".format(
                    element, ion, benchmark.name)) if show_homogenised else None

            figures[benchmark.name] = _compare_repeat_spectra(measurements,
                self.release.node_colors, homogenised_measurements, 
                scaled=scaled, abundance_extent=abundance_extent, bins=bins,
                highlight_flagged=highlight_flagged,
                reference_abundance=benchmark.abundances[element][0],
                reference_uncertainty=benchmark.abundances[element][1],
                reference_label=benchmark.name,
                **kwargs)

            if figures[benchmark.name] is None:
                logger.warn("No finite measurments for benchmark {0}".format(
                    benchmark.name))

        return figures


    def benchmark_line_abundances(self, element, ion, benchmark_filename,
        scaled=False, highlight_flagged=True, show_homogenised=True,
        sort_by=None, differential_abundance_extent=(-1, +1), **kwargs):
        """
        Show the difference in line abundances for each benchmark and each line,
        based on the Jofre et al. (2015) abundances.

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

        :param sort_by: [optional]
            Sort the benchmarks by some parameter (e.g., teff, logg, feh). The
            default is name.

        :type sort_by:
            str

        :param differential_abundance_extent: [optional]
            The lower and upper range to show in differential abundances.

        :type differential_abundance_extent:
            None or a two-length tuple of floats
        
        :returns:
            A figure showing the difference in the line abundances for all of
            the benchmark stars.
        """

        benchmarks = utils.load_benchmarks(benchmark_filename)
        if len(benchmarks) == 0:
            logger.warn("No benchmarks found in {0}".format(benchmark_filename))
            return None

        ok = lambda bm: np.isfinite(bm.abundances.get(element, [np.nan])[0])
        use_benchmarks = filter(ok, benchmarks)

        if len(use_benchmarks) == 0:
            logger.warn("No useable benchmarks. The sample of benchmarks ({0})"\
                " did not contain any stars with measured {1} abundances."\
                .format(len(benchmarks), element))
            return None

        if len(use_benchmarks) < len(benchmarks):
            logger.warn("The following benchmark stars were excluded because "\
                "they had no {0} measurements: {1}".format(element,
                    ", ".join([bm.name for bm in \
                        set(benchmarks).difference([use_benchmarks])])))

        # Sort the useable benchmarks.
        sort_by = sort_by or "name"
        benchmarks = sorted(use_benchmarks, key=lambda bm: getattr(bm, sort_by))
        column = "scaled_abundance" if scaled else "abundance"

        data = self.release.retrieve_table("""SELECT * FROM line_abundances l
            JOIN (SELECT cname, ges_fld, object FROM node_results WHERE ges_type
            LIKE '%_BM') n ON (l.element = '{0}' AND l.ion = '{1}'
            AND l.cname = n.cname AND l.{2} <> 'NaN') ORDER BY wavelength ASC"""\
            .format(element, ion, column))

        nodes = sorted(set(data["node"]))
        wavelengths = sorted(set(data["wavelength"]))

        Nx, Ny = len(nodes), len(wavelengths)

        xscale, yscale, escale = (8, 2, 2)
        xspace, yspace = (0.05, 0.1)
        lb, tr = (0.5, 0.2)
        xs = xscale * Nx + xscale * (Nx - 1) * xspace
        ys = yscale * Ny + yscale * (Ny - 1) * yspace

        xdim = lb * escale + xs + tr * escale
        ydim = lb * escale + ys + tr * escale 

        fig, axes = plt.subplots(Ny, Nx, figsize=(xdim, ydim))
        fig.subplots_adjust(
            left=(lb * escale)/xdim,
            right=(tr * escale + xs)/xdim,
            bottom=(lb * escale)/ydim,
            top=(tr * escale + ys)/ydim,
            wspace=xspace, hspace=yspace)

        scatter_kwds = {
            "s": 50,
            "zorder": 10
        }

        data = data.group_by(["wavelength"])
        for i, (ax_group, wavelength, group) \
        in enumerate(zip(axes, wavelengths, data.groups)):

            x_data = { node: [] for node in nodes }
            y_data = { node: [] for node in nodes }
            is_flagged = { node: [] for node in nodes }
            y_err = { node: [] for node in nodes }
            y_mean_offsets = { node: [] for node in nodes }

            for j, benchmark in enumerate(benchmarks):
                logger.debug("Matching on {}".format(benchmark))

                # Either matches on OBJECT or GES_FLD
                name = benchmark.name.lower()
                group["ges_fld"] = map(str.strip, map(str.lower, group["ges_fld"]))
                group["1.object"] = map(str.strip, map(str.lower, group["1.object"]))
                match = np.array([k for k, row in enumerate(group) \
                    if name == row["ges_fld"] or name == row["1.object"]])

                logger.debug("Found {0} matches for {1}".format(
                    len(match), benchmark))

                if not any(match):
                    logger.warn("Found no benchmark matches for {0}"\
                        .format(benchmark))
                    continue

                for node in nodes:
                    match_node = match[group["node"][match] == node]

                    x_data[node].extend(j * np.ones(len(match_node)))
                    difference = group["abundance"][match_node] \
                        - benchmark.abundances[element][0]
                    y_data[node].extend(difference)
                    y_err[node].extend(group["e_abundance"][match_node])

                    y_mean_offsets[node].append(
                        np.nanmean(difference[group["flags"][match_node] == 0]))
                    is_flagged[node].extend(group["flags"][match_node] > 0)

            # Arrayify!
            for node in nodes:
                x_data[node] = np.array(x_data[node])
                y_data[node] = np.array(y_data[node])
                y_err[node] = np.array(y_err[node])
                is_flagged[node] = np.array(is_flagged[node])

            for k, (ax, node) in enumerate(zip(ax_group, nodes)):
                if ax.is_first_row():
                    ax.set_title(node)

                if ax.is_first_col():
                    ax.text(0.05, 0.95, latexify(wavelength),
                        transform=ax.transAxes, horizontalalignment="left",
                        verticalalignment="top")

                ax.axhline(0, c="k", zorder=0)

                if ax.is_first_col():
                    ax.set_ylabel(r"$\Delta\log_{\epsilon}({\rm X})$")
                else:
                    ax.set_yticklabels([])

                ax.set_xlim(0, len(benchmarks))
                ax.set_xticks(0.5 + np.arange(len(benchmarks)))
                if ax.is_last_row():
                    ax.set_xticklabels([bm.name for bm in benchmarks],
                        rotation=90)
                else:
                    ax.set_xticklabels([])

                color = self.release.node_colors[node]

                flagged = is_flagged[node]
                if highlight_flagged and any(flagged):
                    x = 0.5 + x_data[node][flagged]
                    y = y_data[node][flagged]
                    yerr = y_err[node][flagged]

                    ax.scatter(x, y, facecolor=color, edgecolor="r", lw=2)
                    ax.errorbar(x, y, yerr=yerr, lc="k", ecolor="k", aa=True, 
                        fmt=None, mec="k", mfc="w", ms=6, zorder=1)
            
                else:
                    x = 0.5 + np.array(x_data[Fnode])[~flagged]
                    y = np.array(y_data[node])[~flagged]
                    yerr = np.array(y_err[node])[~flagged]

                    ax.scatter(x, y, facecolor=color, **scatter_kwds)
                    ax.errorbar(x, y, yerr=yerr, lc="k", ecolor="k", aa=True, 
                        fmt=None, mec="k", mfc="w", ms=6, zorder=1)               
                
                # Show relative mean and std. dev for each node
                mean = np.nanmean(y_data[node][flagged])
                sigma = np.nanstd(y_data[node][flagged])

                ax.axhline(np.nanmean(y_mean_offsets[node]), c=color, lw=2,
                    linestyle=":")
                ax.axhspan(mean - sigma, mean + sigma, ec=None, fc=color, 
                    alpha=0.5, zorder=-10)
                ax.axhline(mean, c=color, lw=2, zorder=-1,
                    label=latexify(node.strip()))

        # Common y-axis limits.
        if differential_abundance_extent is None:
            y_lim = max([np.abs(ax.get_ylim()).max() for ax in axes.flatten()])
            [ax.set_ylim(-y_lim, +y_lim) for ax in axes.flatten()]
        else:
            #[ax.set_ylim(extent) for ax in axes.flatten()]
            plt.draw()

            lower, upper = differential_abundance_extent
            # Mark how many measurements are out of frame.
            for ax in axes.flatten():
                data = ax.collections[-1].get_offsets()

                too_high = Counter(data[:, 0][data[:, 1] > upper])
                too_low = Counter(data[:, 0][data[:, 1] < lower])

                print(differential_abundance_extent)
                print(too_high)
                print(too_low)
                if len(too_high) > 0 or len(too_low) > 0:
                    raise a

            raise a

        #  TODO Mark how many lines are out of the frame.

        raise a

        fig.tight_layout()
        return fig



    def transition_heatmap(self, element, ion, column="abundance", **kwargs):
        """
        Display a heatmap of lines that were used for a given species.

        :param element:
            The name of the element to display.

        :type element:
            str

        :param ion:
            The ionisation state of the element (1 = neutral).

        :type ion:
            int

        :param column: [optional]
            The column name to display the heat map for. Default is abundance.

        :type column:
            str

        :param linear: [optional]
            Show a linear scaling in the heatmap. Default is to use logarithmic.

        :type linear:
            bool

        :returns:
            A heatmap showing the usage of individual lines from each node.
        """

        column = column.lower()
        _ = self.release.retrieve_table("SELECT * FROM line_abundances LIMIT 1")
        if column not in _.dtype.names:
            raise ValueError("column '{0}' does not exist".format(column))
        
        try:
            float(_[column])
        except:
            raise ValueError("column '{0} does not represent a float".format(
                column))

        data = self.release.retrieve_table("""SELECT * FROM line_abundances
            WHERE element = %s AND ion = %s AND {0} <> 'NaN'""".format(column),
            (element, ion))

        nodes = sorted(set(data["node"]))
        wavelengths = sorted(set(data["wavelength"]))
        N_nodes, N_wavelengths = map(len, (nodes, wavelengths))

        count = np.zeros((N_nodes, N_wavelengths))
        for i, node in enumerate(nodes):
            match = (data["node"] == node)
            for j, wavelength in enumerate(wavelengths):
                count[i, j] = (match * (data["wavelength"] == wavelength)).sum()

        kwds = {
            "aspect": "auto",
            "cmap": colormaps.plasma,
            "interpolation": "nearest"
        }
        kwds.update(kwargs)

        cube = 0.5
        xs, ys = cube * N_wavelengths, cube * N_nodes
        xspace, yspace = 0, 0
        left, bottom = 1.5, 1.5,
        top, right = 0.2, 0.2
        
        xdim = left + xs + right 
        ydim = bottom + ys + top

        fig, ax = plt.subplots(figsize=(xdim, ydim))
        fig.subplots_adjust(
            left=(left)/xdim,
            right=(left + xs)/xdim,
            bottom=(bottom)/ydim,
            top=(bottom + ys)/ydim,
            wspace=xspace, hspace=yspace)
        image = ax.imshow(count, **kwds)

        ax.set_xlabel(r"Wavelength $[\AA]$")
        
        ax.set_xticks(np.arange(N_wavelengths) + 0.5)
        ax.set_xticklabels(["{0:.1f}".format(_) for _ in wavelengths],
            rotation=90)

        ax.set_yticks(np.arange(N_nodes))
        ax.set_yticklabels(nodes)

        ax.xaxis.set_tick_params(width=0)
        ax.yaxis.set_tick_params(width=0)

        cbar = plt.colorbar(image, ax=[ax])
        cbar.locator = MaxNLocator(4)
        cbar.update_ticks()
        cbar.ax.set_aspect(8)
        cbar.set_label(r"$N$")
        box = cbar.ax.get_position()

        # Adjust the height.        
        required_aspect = float(N_wavelengths)/N_nodes
        box = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        width, height = box.width, box.height

        new_height = width/required_aspect + bottom + top
        fig.set_size_inches(fig.get_figwidth(), new_height, forward=True)

        return fig


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
    if not np.any(np.isfinite(measurements[column])): return None

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

    if show_legend:
        for node in nodes:
            if np.any(np.isfinite(measurements[column]) \
                * (measurements["node"] == node)):
                axes[0][0].plot([], [], color=colors[node], label=node)
        axes[0][0].legend(loc="upper left", frameon=False, fontsize=12)

    # Set S/N axes on the same x-scale.
    xlims = np.array([ax.get_xlim() for ax in axes[:, 1]])
    [ax.set_xlim(xlims[:, 0].min(), xlims[:,1].max()) for ax in axes[:, 1]]

    return fig
