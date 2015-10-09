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
            """SELECT * FROM line_abundances l JOIN
            (SELECT DISTINCT ON (cname) cname, snr FROM node_results
            ORDER BY cname) n ON (trim(l.element) = '{0}'
                AND l.ion = '{1}' AND l.cname = n.cname
                AND l.cname LIKE 'ssssss%')""".format(element, ion))
        if measurements is None: return None

        homogenised_measurements = self.release.retrieve_table(
            """SELECT * FROM homogenised_line_abundances l JOIN
            (SELECT DISTINCT ON (cname) cname, snr 
                FROM node_results ORDER BY cname) n 
            ON (trim(l.element) = '{0}' AND l.ion = '{1}' AND l.cname = n.cname
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
            (trim(l.element) = '{0}' AND l.ion = '{1}' AND l.cname = n.cname 
                AND n.ges_fld LIKE 'M67%' AND n.object = '1194')"""\
            .format(element, ion))
        if measurements is None: return None

        homogenised_measurements = self.release.retrieve_table(
            """SELECT * FROM homogenised_line_abundances l JOIN (SELECT
            DISTINCT ON (cname) cname, snr, ges_fld, object FROM node_results
            ORDER BY cname) n ON (trim(l.element) = '{0}' AND l.ion = '{1}' AND 
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
            WHERE trim(element) = %s AND ion = %s {flag_query}
            ORDER BY wavelength ASC
            """.format(flag_query="AND flags = 0" if ignore_flagged else ""),
            (element, ion))
        if data is None: return

        # Use only finite measurements.
        column = "scaled_abundance" if scaled else "abundance"
        use = np.isfinite(data[column]) * (data["upper_abundance"] == 0)
        if not any(use): return None

        data = data[use]
        wavelengths = sorted(set(data["wavelength"]))

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
        axes = np.atleast_2d(axes)

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
        if X_diff is None:
            return fig

        X_diff = X_diff[np.isfinite(X_diff)]

        b_min, b_max = differential_abundance_extent \
            or (np.nanmin(X_diff), np.nanmax(X_diff))
        hist_kwds["bins"] = np.linspace(b_min, b_max, bins + 1)

        for i, (ax, wavelength) in enumerate(zip(axes.T[1], wavelengths)):
            ax.set_title(latexify(wavelength))
            ax.xaxis.set_major_locator(MaxNLocator(5))
            ax.yaxis.set_major_locator(MaxNLocator(5))
            if not ax.is_last_row():
                ax.set_xticklabels([])
            else:
                ax.set_xlabel(r"$\Delta${0} {1}".format(element, ion))

            match = diff_data["wavelength"] == wavelength
            X_diff_wavelength = utils.calculate_differential_abundances(
                X[match], full_output=False).flatten()

            ax.text(0.95, 0.95, latexify(X_diff.size), transform=ax.transAxes,
                verticalalignment="top", horizontalalignment="right",
                color=full_distribution_color)
            ax.text(0.95, 0.95, 
                latexify("\n{}".format(np.isfinite(X_diff_wavelength).sum())),
                verticalalignment="top", horizontalalignment="right",
                color=comp_distribution_color, transform=ax.transAxes)

            if X_diff.size == 0:
                continue

            if ax.is_first_row():
                ax.text(0.05, 0.95,
                    r"$\mu = {0:.2f}$" "\n" r"$\sigma = {1:.2f}$".format(
                    np.nanmean(X_diff), np.nanstd(X_diff)),  fontsize=10,
                    transform=ax.transAxes, color=full_distribution_color,
                    verticalalignment="top", horizontalalignment="left")
                    
            # Show the full distribution of differential abundances.
            ax.hist(X_diff, color=full_distribution_color, **hist_kwds)

            # Show the distribution of differential abundances for this 
            # wavelength.
            if np.isfinite(X_diff_wavelength).sum() > 0:
                ax.hist(X_diff_wavelength, color=comp_distribution_color,
                    **hist_kwds)
            
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
                
                try:
                    ax.hist(
                        X_diff_wavelength_node[np.isfinite(X_diff_wavelength_node)],
                        color=colors[node], **hist_kwds)
                except:
                    logger.exception("Could not plot differential distribution"\
                        " for {0}".format(node))
                    None

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
                figures[benchmark.name] = None
                logger.warn("No elemental abundance for {0} in {1}".format(
                    element, benchmark.name))
                continue

            measurements = self.release.retrieve_table(
                """SELECT * FROM line_abundances l JOIN 
                (SELECT DISTINCT ON (cname) cname, snr, ges_fld, ges_type, object
                    FROM node_results ORDER BY cname) n
                ON (trim(l.element) = '{0}' AND l.ion = '{1}' AND
                    l.cname = n.cname AND ges_type LIKE '%_BM' AND
                    (ges_fld ILIKE '{2}%' OR n.object ILIKE '{2}%'))""".format(
                    element, ion, benchmark.name))
            if measurements is None:
                figures[benchmark.name] = None
                logger.warn("No {0} {1} measurements for benchmarks".format(
                    element, ion))
                continue

            homogenised_measurements = self.release.retrieve_table(
                """SELECT * FROM homogenised_line_abundances l 
                JOIN (SELECT DISTINCT ON (cname) cname, snr, ges_fld, ges_type,
                    object FROM node_results ORDER BY cname) n
                ON (trim(l.element) = '{0}' AND l.ion = '{1}' AND
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

        data = self.release.retrieve_table(
            """SELECT distinct on (l.abundance_filename, l.wavelength, l.element,
            l.ion) * FROM line_abundances l JOIN (SELECT cname, ges_fld, object
            FROM node_results WHERE ges_type LIKE '%_BM') n 
            ON (trim(l.element) = '{0}' AND l.ion = '{1}' AND l.cname = n.cname 
            AND l.{2} <> 'NaN') ORDER BY l.abundance_filename, l.wavelength,
            l.element, l.ion ASC""".format(element, ion, column))

        homogenised_data = self.release.retrieve_table(
            """SELECT DISTINCT ON (l.spectrum_filename_stub) * FROM
            homogenised_line_abundances l JOIN (SELECT cname, ges_fld,
            object FROM node_results WHERE ges_type like '%_BM') n ON
            (trim(l.element) = '{0}' AND l.ion = '{1}' AND l.cname = n.cname
            AND l.abundance <> 'NaN') ORDER BY l.spectrum_filename_stub"""\
            .format(element, ion)) if show_homogenised else None

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
        axes = np.atleast_2d(axes)

        scatter_kwds = {
            "s": 50,
            "zorder": 10
        }

        homogenised_scatter_kwds = {
            "facecolor": "w",
            "zorder": 100,
            "s": 50,
            "lw": 2
        }

        data = data.group_by(["wavelength"])

        if homogenised_data is not None:
            
            homogenised_x_data = { w: [] for w in wavelengths }
            homogenised_y_data = { w: [] for w in wavelengths }
            homogenised_y_err = { w: [] for w in wavelengths }

            homogenised_data = homogenised_data.group_by(["wavelength"])
            homogenised_data["object"] = \
                map(str.strip, map(str.lower, homogenised_data["object"]))
            homogenised_data["ges_fld"] = \
                map(str.strip, map(str.lower, homogenised_data["ges_fld"]))

            for group in homogenised_data.groups:
                wavelength = group["wavelength"][0]
                for j, benchmark in enumerate(benchmarks):
                    logger.debug("Matching homogenised for {}".format(benchmark))

                    name = benchmark.name.lower()
                    match = np.array([k for k, row in enumerate(group) \
                        if name == row["ges_fld"] or name == row["object"]])
                    logger.debug("Found {0} homogenised matches for {1}".format(
                        len(match), benchmark.name))

                    if len(match) == 0: continue

                    homogenised_x_data[wavelength].extend(
                        j * np.ones(len(match)))
                    homogenised_y_data[wavelength].extend(
                        group["abundance"][match] \
                            - benchmark.abundances[element][0])
                    homogenised_y_err[wavelength].extend(
                        group["e_abundance"][match])
                
                # Arrayify!
                
            for wavelength in wavelengths:
                homogenised_x_data[wavelength] \
                    = np.array(homogenised_x_data[wavelength])
                homogenised_y_data[wavelength] \
                    = np.array(homogenised_y_data[wavelength])
                homogenised_y_err[wavelength] \
                    = np.array(homogenised_y_err[wavelength])

        else:
            homogenised_x_data = {}
            homogenised_y_data = {}
            homogenised_y_err = {}

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
                # TODO do this better
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
                    difference = group[column][match_node] \
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
                ax.set_xticks(np.arange(len(benchmarks)))
                if ax.is_last_row():
                    ax.set_xticklabels([bm.name for bm in benchmarks],
                        rotation=90, verticalalignment="top",
                        horizontalalignment="left")
                else:
                    ax.set_xticklabels([])

                if len(x_data[node]) == 0: continue

                color = self.release.node_colors[node]

                flagged = is_flagged[node]
                if highlight_flagged and any(flagged):
                    x = 0.5 + x_data[node][flagged]
                    y = y_data[node][flagged]
                    yerr = y_err[node][flagged]

                    ax.scatter(x, y, facecolor=color, edgecolor="r", lw=2)
                    ax.errorbar(x, y, yerr=yerr, lc="k", ecolor="k", aa=True, 
                        fmt=None, mec="k", mfc="w", ms=6, zorder=1)
            
                x = 0.5 + np.array(x_data[node])[~flagged]
                y = np.array(y_data[node])[~flagged]
                yerr = np.array(y_err[node])[~flagged]

                ax.scatter(x, y, facecolor=color, **scatter_kwds)
                ax.errorbar(x, y, yerr=yerr, lc="k", ecolor="k", aa=True, 
                    fmt=None, mec="k", mfc="w", ms=6, zorder=9)               

                # Show relative mean and std. dev for each node
                mean = np.nanmean(y_data[node][~flagged])
                sigma = np.nanstd(y_data[node][~flagged])

                ax.axhline(np.nanmean(y_mean_offsets[node]), c=color, lw=2,
                    linestyle=":")
                ax.axhspan(mean - sigma, mean + sigma, ec=None, fc=color, 
                    alpha=0.5, zorder=-10)
                ax.axhline(mean, c=color, lw=2, zorder=-1,
                    label=latexify(node.strip()))

            # Draw all homogenised values, if they exist.
            if homogenised_data is not None:
                x = 0.5 + homogenised_x_data[wavelength]
                y = homogenised_y_data[wavelength]
                yerr = homogenised_y_err[wavelength]
                
                for ax in ax_group:
                    ax.errorbar(x, y, yerr=yerr, lc="k", ecolor="k", aa=True,
                        fmt=None, mec="k", mfc="w", ms=6, lw=2, zorder=99)
                    ax.scatter(x, y, **homogenised_scatter_kwds)

            if differential_abundance_extent is not None:
                extent = lower, upper = differential_abundance_extent
                
                y_high = 0.95 * np.ptp(extent) + lower
                y_low = 0.05 * np.ptp(extent) + lower

                for ax, node in zip(ax_group, nodes):

                    ax.set_ylim(extent)

                    too_high = Counter(x_data[node][y_data[node] > upper])
                    too_low = Counter(x_data[node][y_data[node] < lower])

                    for x_position, N in too_high.items():
                        ax.text(x_position + .5, y_high, latexify(N), color="r",
                            fontsize=10, zorder=100, verticalalignment="top",
                            horizontalalignment="center")

                    for x_position, N in too_low.items():
                        ax.text(x_position + .5, y_low, latexify(N), color="r",
                            fontsize=10, zorder=100, verticalalignment="bottom",
                            horizontalalignment="center")

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
            WHERE trim(element) = %s
            AND ion = %s AND {0} <> 'NaN'""".format(column), (element, ion))
        if data is None: return

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
        
        ax.set_xticks(np.arange(N_wavelengths))
        ax.set_xticklabels(["{0:.1f}".format(_) for _ in wavelengths],
            rotation=90)

        ax.set_yticks(range(N_nodes))
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


    def mean_abundances(self, element, ion, bins=None, extent=None, **kwargs):
        """
        Show the mean (node-reported) abundances of the given species from each
        node.

        :param element:
            The element to show the mean abundances for.

        :type element:
            str

        :param ion:
            The ionisation stage of the element to show mean abundances for (an
            ionisation stage of 1 indicates neutral).

        :type ion:
            int

        :param bins: [optional]
            The bins to use in the histograms.

        :type bins:
            int or array of bin edges

        :param extent: [optional]
            The range in abundances to show.

        :type extent:
            None or two-length tuple

        :returns:
            A figure showing the mean abundances (as reported by the nodes), as
            compared to all other nodes.
        """

        # Check that the column is actually in the node results table!
        column = "{0}{1}".format(element.lower(), ion)
        check = self.release.retrieve_table("SELECT * FROM node_results LIMIT 1")
        if column not in check.dtype.names:
            logger.warn("Column {0} not in node_results table".format(column))
            return None

        nodes = self.release.nodes
        N_nodes = len(nodes)
        N_entries = int(self.release.retrieve("""SELECT count(*) FROM
            node_results WHERE node = %s""", (nodes[0], ))[0][0])

        X = np.nan * np.ones((N_nodes, N_entries))
        X_uncertainties = X.copy()
        
        
        for i, node in enumerate(nodes):
            data = self.release.retrieve_table("""SELECT {0}, e_{0}
                FROM node_results WHERE node = %s ORDER BY ra DESC, dec DESC,
                cname DESC""".format(column), (node, ))

            if data is None:
                logger.warn("No node results found for node {0}".format(node))
                continue
            
            X[i, :] = data[column]
            X_uncertainties[i, :] = data["e_{}".format(column)]
            
        if 2 > np.max(np.sum(np.isfinite(X), axis=0)):
            # No data for any star from 2 nodes to compare against.
            return None

        bins = bins or np.arange(-0.50, 0.55, 0.05)
        return _corner_scatter(X, uncertainties=X_uncertainties,
            labels=map(str.strip, nodes), bins=bins, extent=extent, **kwargs)


    def line_abundances(self, element, ion, reference_column, scaled=False,
        highlight_flagged=True, show_homogenised=True, aux_column=None,
        x_extent=None, y_extent=None, show_node_comparison=True,
        show_line_comparison=True, show_legend=False, show_uncertainties=False,
        abundance_format="log_x", **kwargs):
        """
        Show the reference and relative abundances for a given element and ion
        against the reference column provided.
        
        :param element:
            The atomic element of interest.

        :type element:
            str

        :param ion:
            The ionisation stage of the species of interest (1 indicates neutral)

        :type ion:
            int

        :param reference_column:
            The name of the reference column (from the node results or line
            abundances tables) to display on the x-axis.

        :type reference_column:
            str

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

        :param aux_column: [optional]
            Show an auxiliary column as a colorbar.

        :type aux_column:
            str

        :param x_extent: [optional]
            The range of values to display in the x-axis.

        :type x_extent:
            None or two-length tuple

        :param y_extent: [optional]
            The range of values to display in the y-axis.

        :type y_extent:
            None or two-length tuple

        :param show_node_comparison: [optional]
            Show triangle markers indicating all abundances from the given node.

        :type show_node_comparison:
            bool

        :param show_line_comparison: [optional]
            Show square markers indicating all abundances for the given line.

        :type show_line_comparison:
            bool

        :param show_legend: [optional]
            Show a legend for the node and line comparisons. This will only be
            displayed if both the show_line_comparison and show_node_comparison
            options are also set to True.

        :type show_legend:
            bool

        :param abundance_format: [optional]
            The abundance format to display. Available options are 'log_x',
            'x_h', or 'x_fe'.

        :type abundance_format:
            str
        """

        # Ensure the reference column is valid.
        reference_column = reference_column.lower()
        _ = self.release.retrieve_table("SELECT * FROM node_results LIMIT 1")
        if reference_column not in _.dtype.names:
            raise ValueError(
                "reference column '{0}' not valid (acceptable: {1})".format(
                    reference_column, ", ".join(check.dtype.names)))

        # Check the abundance format.
        available = ("x_h", "x_fe", "log_x")
        abundance_format = abundance_format.lower()
        if abundance_format not in available:
            raise ValueError(
                "abundance format '{0}' is not valid (acceptable: {1})".format(
                    abundance_format, ", ".join(available)))

        columns = ["feh", reference_column]
        if aux_column is not None: columns += [aux_column]
        
        data = self.release.retrieve_table(
            """SELECT * FROM line_abundances l JOIN 
            (SELECT DISTINCT ON (cname) cname, {0} FROM node_results 
                ORDER BY cname) n
            ON (trim(l.element) = %s AND l.ion = %s AND l.cname = n.cname)
            ORDER BY abundance_filename, wavelength ASC""".format(
                ", ".join(set(columns))), (element, ion))
        if data is None: return            

        scatter_kwds = {
            "s": 50
        }
        if aux_column is not None:
            scatter_kwds["cmap"] = colormaps.plasma
            aux_extent = kwargs.pop("aux_extent", None)
            if aux_extent is None:
                scatter_kwds["vmin"] = np.nanmin(data[aux_column])
                scatter_kwds["vmax"] = np.nanmax(data[aux_column])
            else:
                scatter_kwds["vmin"], scatter_kwds["vmax"] = aux_extent

        homogenised_scatter_kwds = scatter_kwds.copy()
        homogenised_scatter_kwds["lw"] = 2

        nodes = sorted(set(data["node"]))
        if show_homogenised:
            assert "Homogenised" not in nodes, "Restricted name"
            homogenised_data = self.release.retrieve_table(
                """SELECT * FROM homogenised_line_abundances l JOIN (SELECT
                DISTINCT ON (cname) cname, {0} FROM node_results ORDER BY cname)
                n ON (trim(l.element) = %s AND l.ion = %s AND l.cname = n.cname)
                ORDER BY wavelength ASC""".format(", ".join(set(columns))),
                    (element, ion))
            if homogenised_data is None:
                show_homogenised = False
            else:
                nodes += ["Homogenised"]

        wavelengths = sorted(set(data["wavelength"]))
        N_nodes, N_lines = len(nodes), len(wavelengths)

        # How many figures should we have?
        axes_per_figure = kwargs.pop("axes_per_figure", 20)
        if np.isfinite(axes_per_figure):
            N_figures \
                = int(np.ceil(float(N_lines) / axes_per_figure))
        else:
            N_figures, axes_per_figure = 1, N_lines

        figures = {}
        for i in range(N_figures):
            w = wavelengths[i*axes_per_figure:(i + 1)*axes_per_figure]
            figures[str(i)] = _line_abundances(element, ion, reference_column,
                data, homogenised_data, nodes, w, scatter_kwds=scatter_kwds,
                homogenised_scatter_kwds=homogenised_scatter_kwds,
                scaled=scaled, highlight_flagged=highlight_flagged,
                show_homogenised=show_homogenised,
                show_node_comparison=show_node_comparison,
                show_line_comparison=show_line_comparison,
                show_legend=show_legend, show_uncertainties=show_uncertainties,
                abundance_format=abundance_format, **kwargs)

        return figures if N_figures > 1 else figures.values()[0]



    def differential_line_abundances_wrt_parameter(self, element, ion, parameter,
        scaled=False, ignore_flagged=False, logarithmic=True, bins=25,
        x_extent=None, y_extent=(-0.5, 0.5), **kwargs):
        """
        Show the node differential line abundances for a given element and ion
        with respect to a given column in the node results or line abundances
        tables.

        :param element:
            The elemental abundance of interest.

        :type element:
            str

        :param ion:
            The ionisation stage of the element of interest (1 = neutral).

        :type ion:
            int

        :param parameter:
            The x-axis parameter to display differential abundances against.

        :type parameter:
            str

        :param scaled: [optional]
            Show scaled abundances. Alternative is to show unscaled abundances.

        :type scaled:
            bool

        :param highlight_flagged: [optional]
            Show a red outline around measurements that are flagged.

        :type highlight_flagged:
            bool

        :param logarithmic: [optional]
            Display logarithmic counts.

        :type logarithmic:
            bool

        :param bins: [optional]
            The number of bins to have in the differential abundances axes. The
            default bin number is 50.

        :type bins:
            int

        :param x_extent: [optional]
            The lower and upper range of the x-axis values to display.

        :type x_extent:
            None or two-length tuple of floats

        :param y_extent: [optional]
            The lower and upper range in differential abundances to display.

        :type y_extent:
            two-length tuple of floats

        :returns:
            A figure showing the differential line abundances against the 
            requested parameter.
        """

        # Get the full distribution of abundances.
        X, nodes, data_table = self.release._match_species_abundances(element,
            ion, additional_columns=[parameter], scaled=scaled,
            include_flagged_lines=not ignore_flagged)
        if X is None: return None

        if len(nodes) == 1:
            logger.warn(
                "Only one node present; differential abundances are unavailable")
            return None

        # Calculate differential abundances.
        differential_abundances, indices \
            = utils.calculate_differential_abundances(X, full_output=True)
        if differential_abundances is None: return None

        wavelengths = [w for w in sorted(set(data_table["wavelength"])) if np.any( \
            np.isfinite(differential_abundances[data_table["wavelength"] == w]))]
        N_nodes, N_wavelengths = len(nodes), len(wavelengths)

        N_original_wavelengths = len(set(data_table["wavelength"]))
        if N_original_wavelengths > N_wavelengths:
            logger.warn("{0} of {1} wavelengths not shown because there are no"\
                " differential measurements".format(
                    N_original_wavelengths - N_wavelengths,
                    N_original_wavelengths))

        if 0 in (N_nodes, N_wavelengths):
            return None

        # Create a figure of the right dimensions.
        Nx, Ny = 1 + N_nodes, N_wavelengths
        xscale, yscale, escale = (4, 2, 2)
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
        axes = np.atleast_2d(axes)

        # Get the common bin sizes.
        x_min, x_max = x_extent or (
            np.floor(np.nanmin(data_table[parameter].astype(float))),
            np.ceil(np.nanmax(data_table[parameter].astype(float)))
        )
        x_bins = np.linspace(x_min, x_max, bins + 1)

        differential_min, differential_max = y_extent or \
            (np.nanmin(differential_abundances), np.nanmax(differential_abundances))
        y_bins = np.linspace(differential_min, differential_max, bins + 1)

        # Get common boundaries for the x- and y-axis.
        histogram_kwds = {
            "normed": kwargs.pop("normed", False),
            "bins": (x_bins, y_bins)
        }
        imshow_kwds = {
            "interpolation": "nearest",
            "aspect": "auto",
            "cmap": colormaps.plasma,
            "extent": (x_bins[0], x_bins[-1], y_bins[0], y_bins[-1]),
        }
        imshow_kwds.update(kwargs)
        axes = np.atleast_2d(axes)

        for i, (row_axes, wavelength) in enumerate(zip(axes, wavelengths)):

            # First one should contain *everything* in a greyscale.
            full_ax, node_axes = row_axes[0], row_axes[1:]

            # For the full axes we have to repeat the x-axis data.
            mask = data_table["wavelength"] == wavelength
            x = np.tile(data_table[parameter].astype(float)[mask],
                differential_abundances.shape[1])
            y = differential_abundances[mask, :].T.flatten()

            H, xe, ye = np.histogram2d(x, y, **histogram_kwds)
            Z = np.log(1 + H.T) if logarithmic else H.T
            full_ax.imshow(Z, **imshow_kwds)
            
            # Bells and whistles
            _ = "{0}\,{1}".format(element, ion)
            full_ax.set_ylabel(
                r"$\Delta\log_{\epsilon}({\rm " + _ + "})$")
            if full_ax.is_first_row():
                full_ax.set_title("All nodes")
            if full_ax.is_last_row():
                full_ax.set_xlabel(parameter)
            else:
                full_ax.set_xticklabels([])

            full_ax.text(0.05, 0.95, r"${0}$ $\AA$".format(wavelength),
                transform=full_ax.transAxes, color="w", fontsize=14,
                verticalalignment="top", horizontalalignment="left")

            # Make the node-specific histogram plots.
            for j, (ax, node) in enumerate(zip(node_axes, nodes)):
                # Some formatting:
                ax.set_yticklabels([])
                if ax.is_first_row():
                    ax.set_title(node.strip())
                if ax.is_last_row():
                    ax.set_xlabel(parameter)
                else:
                    ax.set_xticklabels([])

                # Get all the differential abundances for this node.
                # Recall Nx = 1 + len(nodes) t.f. len(nodes) - 1 = Nx - 2
                node_x \
                    = np.tile(data_table[parameter].astype(float)[mask], Nx - 2)
                if len(node_x) == 0:
                    continue
                
                # Need the node_y abundances to be log_x(NODE A) - log_x(ALL)
                node_y = np.hstack([
                    [-1, +1][j == idx[0]] * differential_abundances[mask, k] \
                    for k, idx in enumerate(indices) if j in idx])

                H, xe, ye = np.histogram2d(node_x, node_y, **histogram_kwds)
                Z = np.log(1 + H.T) if logarithmic else H.T
                image = ax.imshow(Z, **imshow_kwds)

            # Add a zero-differential marker for all lines.
            [ax.axhline(0, c="k", lw=1) for ax in row_axes]

        for ax in axes.flatten():
            ax.xaxis.set_major_locator(MaxNLocator(5))
            ax.yaxis.set_major_locator(MaxNLocator(5))
        
        return fig


    def homogenised_line_abundance_uncertainties(self, element, ion, bins=50,
        x_extent=None):
        """
        Draw histograms showing the distribution of uncertainties in each line
        abundance for the given element.

        :param element:
            The elemental abundance of interest.

        :type element:
            str

        :param ion:
            The ionisation stage of the element of interest (1 = neutral).

        :type ion:
            int

        :param bins: [optional]
            The number of bins to have. The default number is 20.

        :type bins:
            int

        :param x_extent: [optional]
            The range of values to display in the x-axis.

        :type x_extent:
            None or two-length tuple

        :returns:
            A single-panel figure showing the distributions of homogenised
            abundance uncertainties for the given element.
        """

        # Get all of the homogenised line data for this element and ion.
        data = self.release.retrieve_table(
            """SELECT wavelength, e_abundance FROM homogenised_line_abundances
            WHERE TRIM(element) = %s AND ion = %s AND e_abundance <> 'NaN'
            AND upper_abundance = 0 ORDER BY wavelength ASC""", (element, ion))
        if data is None: return None

        data = data.group_by(["wavelength"])

        extent = x_extent \
            or (np.nanmin(data["e_abundance"]), np.nanmax(data["e_abundance"]))
        bins = np.linspace(extent[0], extent[1], 1 + bins)

        hist_kwargs = {
            "color": "#666666",
            "histtype": "step",
            "bins": bins,
            "normed": True
        }
        fig, ax = plt.subplots(1)

        for group in data.groups:
            if np.any((group["e_abundance"] > hist_kwargs["bins"][0]) * \
                (group["e_abundance"] < hist_kwargs["bins"][-1])):
                ax.hist(group["e_abundance"], **hist_kwargs)
            else:
                logger.warn("Skipping over {0} {1} elements because they are "\
                    "outside the display range".format(element, ion))

        # Show the full distribution as a thick line.
        hist_kwargs.update({
            "lw": 3,
            "color": "k"
        })
        hist_kwargs["lw"] = 3

        if np.any((data["e_abundance"] > hist_kwargs["bins"][0]) * \
            (data["e_abundance"] < hist_kwargs["bins"][-1])):
            ax.hist(data["e_abundance"], **hist_kwargs)

        else:
            logger.warn("ALL ELEMENTAL ABUNDANCE UNCERTAINTIES FOR {0} {1} ARE"\
                " OUTSIDE THE DISPLAY RANGE".format(element, ion))

        ax.xaxis.set_major_locator(MaxNLocator(5))
        ax.yaxis.set_major_locator(MaxNLocator(5))

        ax.set_xlabel(r"$\sigma_{" + element + "\,{0}".format(ion) \
            + r"}$ $({\rm dex})$")
        ax.set_ylabel("Normalised count")

        return fig


    def homogenised_abundance_uncertainties(self, element, ion, bins=50,
        x_extent=None):
        """
        Draw histograms showing the distribution of uncertainties in each star's
        abundance for the given element.

        :param element:
            The elemental abundance of interest.

        :type element:
            str

        :param ion:
            The ionisation stage of the element of interest (1 = neutral).

        :type ion:
            int

        :param bins: [optional]
            The number of bins to have. The default number is 20.

        :type bins:
            int

        :param x_extent: [optional]
            The range of values to display in the x-axis.

        :type x_extent:
            None or two-length tuple

        :returns:
            A single-panel figure showing the distributions of homogenised
            abundance uncertainties for the given element.
        """

        # Get all of the homogenised line data for this element and ion.
        data = self.release.retrieve_table(
            """SELECT e_abundance FROM homogenised_abundances
            WHERE TRIM(element) = %s AND ion = %s AND e_abundance <> 'NaN'
            AND upper_abundance = 0""",
            (element, ion))
        if data is None: return None

        extent = x_extent \
            or (np.nanmin(data["e_abundance"]), np.nanmax(data["e_abundance"]))
        bins = np.linspace(extent[0], extent[1], 1 + bins)

        hist_kwargs = {
            "color": "#666666",
            "histtype": "step",
            "bins": bins,
            "normed": True
        }
        fig, ax = plt.subplots(1)

        for group in data.groups:
            if np.any((group["e_abundance"] > hist_kwargs["bins"][0]) * \
                (group["e_abundance"] < hist_kwargs["bins"][-1])):
                ax.hist(group["e_abundance"], **hist_kwargs)
            else:
                logger.warn("Skipping over {0} {1} elements because they are "\
                    "outside the display range".format(element, ion))

        # Show the full distribution as a thick line.
        hist_kwargs.update({
            "lw": 3,
            "color": "k"
        })
        hist_kwargs["lw"] = 3

        if np.any((data["e_abundance"] > hist_kwargs["bins"][0]) * \
            (data["e_abundance"] < hist_kwargs["bins"][-1])):
            ax.hist(data["e_abundance"], **hist_kwargs)

        else:
            logger.warn("ALL ELEMENTAL ABUNDANCE UNCERTAINTIES FOR {0} {1} ARE"\
                " OUTSIDE THE DISPLAY RANGE".format(element, ion))

        ax.xaxis.set_major_locator(MaxNLocator(5))
        ax.yaxis.set_major_locator(MaxNLocator(5))

        ax.set_xlabel(r"$\sigma_{" + element + "\,{0}".format(ion) \
            + r"}$ $({\rm dex})$")
        ax.set_ylabel("Normalised count")

        return fig



def latexify(label):
    """ A placeholder for a smart function to latexify common labels. """
    return label


def _corner_scatter(data, labels=None, uncertainties=None, extent=None,
    color=None, bins=20, relevant_text=None):
    """
    Create a corner scatter plot showing the differences between each node.

    :param extent: [optional]
        The (minimum, maximum) extent of the plots.

    :type extent:
        two-length tuple

    :returns:
        A matplotlib figure.
    """

    # How many nodes to plot?
    N = data.shape[0]
    K = N
    assert K > 0, "Need more than one node to compare against."

    factor = 2.0           # size of one side of one panel
    lbdim = 0.5 * factor   # size of left/bottom margin
    trdim = 0.5 * factor   # size of top/right margin
    whspace = 0.15         # w/hspace size
    plotdim = factor * K + factor * (K - 1.) * whspace
    dim = lbdim + plotdim + trdim

    fig, axes = plt.subplots(K, K, figsize=(dim, dim))
    axes = np.atleast_2d(axes)

    lb = lbdim / dim
    tr = (lbdim + plotdim) / dim
    fig.subplots_adjust(left=lb, bottom=lb, right=tr, top=tr,
        wspace=whspace, hspace=whspace)

    hist_kwargs = {
        "color": "k",
        "histtype": "step"
    }

    extent = extent or (0.9 * np.nanmin(data), 1.1 * np.nanmax(data))
    
    # Match all of the nodes
    for i in range(N):
        for j in range(N):

            if j > i: # hide.
                try:
                    ax = axes[i, j]
                    ax.set_frame_on(False)
                    if ax.is_last_col() and ax.is_first_row() and relevant_text:
                        ax.text(0.95, 0.95, relevant_text, fontsize=14,
                            verticalalignment="top", horizontalalignment="right")
                        [_([]) for _ in (ax.set_xticks, ax.set_yticks)]
                    else:
                        ax.set_visible(False)
                    
                except IndexError:
                    None
                continue
                
            ax = axes[i, j]

            if i == j:
                indices = np.arange(N)
                indices = indices[indices != i]

                diff = (data[i] - data[indices]).flatten()
                diff = diff[np.isfinite(diff)]
                if diff.size and any((bins[-1] > diff) * (diff > bins[0])):
                    ax.hist(diff, bins=bins, **hist_kwargs)
                
            else:    
                ax.plot(extent, extent, "k:", zorder=-100)
                if uncertainties is not None:
                    ax.errorbar(data[i], data[j], 
                        xerr=uncertainties[i], yerr=uncertainties[j],
                        ecolor="k", aa=True, fmt=None, mec='k', mfc="w", ms=6,
                        zorder=1, lc="k")

                ax.scatter(data[i], data[j], c=color or "w", zorder=100)
                
                ax.set_xlim(extent)
                ax.set_ylim(extent)

            ax.xaxis.set_major_locator(MaxNLocator(5))
            ax.yaxis.set_major_locator(MaxNLocator(5))

            if i < K - 1:
                ax.set_xticklabels([])
            else:
                [l.set_rotation(45) for l in ax.get_xticklabels()]
                if labels is not None:
                    if i == j:
                        ax.set_xlabel("{} $-$ $X$".format(labels[j]))
                    else:
                        ax.set_xlabel(labels[j])
                    ax.xaxis.set_label_coords(0.5, -0.3)
                
            if j > 0 or i == j:
                ax.set_yticklabels([])
            else:
                [l.set_rotation(45) for l in ax.get_yticklabels()]
                if labels is not None:
                    ax.set_ylabel(labels[i])
                    ax.yaxis.set_label_coords(-0.3, 0.5)

    return fig


def _line_abundances(element, ion, reference_column, 
    data, homogenised_data, nodes, wavelengths,  
    scatter_kwds=None, homogenised_scatter_kwds=None, scaled=False,
    highlight_flagged=True, show_homogenised=True, aux_column=None,
    x_extent=None, y_extent=None, show_node_comparison=True,
    show_line_comparison=True, show_legend=False, show_uncertainties=False,
    abundance_format="log_x", **kwargs):

    # Calculate figure size
    N_nodes, N_lines = len(nodes), len(wavelengths)
    xdim, scale = None, 1.
    while xdim is None or xdim > 400:

        lb, tr = 0.5 * scale, 0.2 * scale
        xscale, yscale = 4.0 * scale, 1.5 * scale
        wspace, hspace = 0.05 * scale, 0.10 * scale
        
        xs = xscale * N_lines + xscale * (N_lines - 1) * wspace
        ys = yscale * N_nodes + yscale * (N_nodes - 1) * hspace
        x_aux = 0 if aux_column is None else 0.5 * xscale + 4 * wspace

        xdim = lb * xscale + xs + x_aux + tr * xscale
        ydim = lb * yscale + ys + tr * yscale

        # In case xdim > 400
        scale -= 0.1

    logger.debug("Requesting figure size: {0}, {1}".format(xdim, ydim))
    fig, axes = plt.subplots(N_nodes, N_lines, figsize=(xdim, ydim))
    fig.subplots_adjust(
        left=(lb * xscale)/xdim,
        bottom=(lb * yscale)/ydim,
        right=(lb * xscale + xs + x_aux)/xdim,
        top=(tr * yscale + ys)/ydim,
        wspace=wspace, hspace=hspace)
    axes = np.atleast_2d(axes)
    if axes.shape[0] == 1: axes = axes.T

    for i, (node_axes, wavelength) in enumerate(zip(axes.T, wavelengths)):

        match_wavelength = (data["wavelength"] == wavelength)
        if show_homogenised:
            match_homogenised_wavelength \
                = (homogenised_data["wavelength"] == wavelength)

        assert len(node_axes) == len(nodes)
        for j, (ax, node) in enumerate(zip(node_axes, nodes)):

            # Labels and titles
            if ax.is_last_row():
                ax.set_xlabel(reference_column)
            if ax.is_first_col():
                ax.set_ylabel(node)               
            if ax.is_first_row():
                ax.set_title(wavelength)

            # Get all measurements for this line by this node.
            if node == "Homogenised":
                ax_data = homogenised_data
                match = match_homogenised_wavelength
                
            else:
                ax_data = data
                match = match_wavelength * (data["node"] == node)

            x = ax_data[reference_column][match]
            y = ax_data["abundance"][match]                 
            if len(x) == 0:
                continue

            if abundance_format == "x_h":
                y -= utils.solar_abundance(element)
            elif abundance_format == "x_fe":
                y -= utils.solar_abundance(element) \
                    + ax_data["feh"][match].astype(float)

            if node == "Homogenised":
                kwds = homogenised_scatter_kwds.copy()
                if aux_column is not None:
                    kwds["c"] = ax_data[aux_column][match]

                scat = ax.scatter(x, y, **kwds)
            else:
                kwds = scatter_kwds.copy()
                    
                if highlight_flagged:
                    flagged = ax_data["flags"][match] > 0
                    if not all(flagged):
                        if aux_column is not None:
                            kwds["c"] = ax_data[aux_column][match][~flagged]
                        scat = ax.scatter(x[~flagged], y[~flagged], **kwds)

                    if any(flagged):
                        if aux_column is not None:
                            kwds["c"] = ax_data[aux_column][match][flagged]
                        kwds["edgecolor"] = "r"
                        kwds["lw"] = 2
                        ax.scatter(x[flagged], y[flagged], **kwds)

                else:
                    scat = ax.scatter(x, y, **kwds)

            if show_uncertainties:
                raise NotImplementedError("sorry, no time")

    comparison_scatter_kwds = {
        "s": 25,
        "c": "#EEEEEE",
        "zorder": -1,
        "marker": "v",
        "alpha": 0.5,
        "edgecolor": "#BBBBBB"
    }
    # We don't need the homogenised node to see comparisons, so let's remove
    # it from the rest.
    if show_homogenised: nodes.remove("Homogenised")

    # Show all other measurements for this line (from all nodes).
    if show_line_comparison: 
        for i, (node_axes, wavelength) in enumerate(zip(axes.T, wavelengths)):
            match = (data["wavelength"] == wavelength)
            x = data[reference_column][match]
            y = data["abundance"][match]
            if abundance_format == "x_h":
                y -= utils.solar_abundance(element)
            elif abundance_format == "x_fe":
                y -= utils.solar_abundance(element) \
                    + data["feh"][match].astype(float)

            if len(x) == 0: continue

            for j, ax in enumerate(node_axes):
                ax.scatter(x, y, label="Common line", 
                    **comparison_scatter_kwds)

    comparison_scatter_kwds["marker"] = "s"
    if show_node_comparison:
        for i, (line_axes, node) in enumerate(zip(axes, nodes)):
            match = (data["node"] == node)
            x = data[reference_column][match]
            y = data["abundance"][match]
            if abundance_format == "x_h":
                y -= utils.solar_abundance(element)
            elif abundance_format == "x_fe":
                y -= utils.solar_abundance(element) \
                    + data["feh"][match].astype(float)

            if len(x) == 0: continue

            for j, ax in enumerate(line_axes):
                ax.scatter(x, y, label="Common node", 
                    **comparison_scatter_kwds)

    if show_node_comparison and show_line_comparison and show_legend:
        axes.T[0][0].legend(loc="upper left", frameon=True, fontsize=12)

    # Common x- and y-axis limits.
    x_limits = [+np.inf, -np.inf]
    y_limits = [+np.inf, -np.inf]
    for ax in axes.flatten():
        if sum([_.get_offsets().size for _ in ax.collections]) == 0:
            continue

        proposed_x_limits = ax.get_xlim()
        proposed_y_limits = ax.get_ylim()

        if proposed_x_limits[0] < x_limits[0]:
            x_limits[0] = proposed_x_limits[0]
        if proposed_y_limits[0] < y_limits[0]:
            y_limits[0] = proposed_y_limits[0]
        if proposed_x_limits[1] > x_limits[1]:
            x_limits[1] = proposed_x_limits[1]
        if proposed_y_limits[1] > y_limits[1]:
            y_limits[1] = proposed_y_limits[1]

    y_limits = y_extent or y_limits
    x_limits = x_extent or x_limits
    for ax in axes.flatten():
        ax.set_xlim(x_limits)
        ax.set_ylim(y_limits)
        ax.xaxis.set_major_locator(MaxNLocator(5))
        ax.yaxis.set_major_locator(MaxNLocator(5))

        if not ax.is_last_row(): ax.set_xticklabels([])
        if not ax.is_first_col(): ax.set_yticklabels([])

    if aux_column is not None:
        cbar = plt.colorbar(scat, ax=list(axes.flatten()))
        cbar.set_label(aux_column)
        _ = axes.T[-1][0].get_position().bounds
        cbar.ax.set_position([
            (lb*xscale + xs + 4*wspace)/xdim, 
            axes.T[-1][-1].get_position().bounds[1],
            (lb*xscale + xs + x_aux)/xdim,
            axes.T[-1][0].get_position().y1 - \
                axes.T[-1][-1].get_position().y0
            ])

        fig.subplots_adjust(
            left=(lb * xscale)/xdim,
            bottom=(lb * yscale)/ydim,
            right=(lb * xscale + xs)/xdim,
            top=(tr * yscale + ys)/ydim,
            wspace=wspace, hspace=hspace)

    return fig


def _compare_repeat_spectra(measurements, colors, homogenised_measurements=None,
    scaled=False, abundance_extent=None, bins=20, highlight_flagged=True,
    show_legend=True, x_column="SNR", reference_abundance=None, 
    reference_uncertainty=None, reference_label=None, **kwargs):

    nodes = sorted(set(measurements["node"]))
    wavelengths = sorted(set(measurements["wavelength"]))
    N_wavelengths = len(wavelengths)

    axes_per_figure = kwargs.pop("axes_per_figure", 50)
    if np.isfinite(axes_per_figure):
        N_figures = int(np.ceil(float(N_wavelengths) / axes_per_figure))
        logger.info("Producing {0} figures with {1} wavelengths per figure"\
            .format(N_figures, axes_per_figure))

    else:
        N_figures, axes_per_figure = 1, N_wavelengths

    figures = {}
    for i in range(N_figures):
        figures[str(i)] = __compare_repeat_spectra(measurements, colors,
            wavelengths[i*axes_per_figure:(i + 1)*axes_per_figure], nodes,
            homogenised_measurements=homogenised_measurements, scaled=scaled,
            abundance_extent=abundance_extent, bins=bins,
            highlight_flagged=highlight_flagged, show_legend=show_legend,
            x_column=x_column, reference_abundance=reference_abundance,
            reference_uncertainty=reference_uncertainty,
            reference_label=reference_label, **kwargs)

    return figures if N_figures > 1 else figures.values()[0]


def __compare_repeat_spectra(measurements, colors, wavelengths, nodes,
    homogenised_measurements=None, scaled=False, abundance_extent=None, bins=20,
    highlight_flagged=True, show_legend=True, x_column="SNR",
    reference_abundance=None, reference_uncertainty=None, reference_label=None,
    **kwargs):

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
    axes = np.atleast_2d(axes)

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
