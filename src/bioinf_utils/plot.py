import math
from collections import Counter

import numpy
import pandas

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib
import matplotlib.patches as mpatches


DARK_GRAY = [0.5] * 3


class _Figure:
    def __init__(self, fig_size=None):
        self._fig = Figure(figsize=fig_size)
        self._canvas = FigureCanvas(
            self._fig
        )  # Don't remove it or savefig will fail later

    def save_fig(self, path):
        self._fig.savefig(path)


class SimpleFigure(_Figure):
    def __init__(self, fig_size=None, subplot_kw: dict = None):
        super().__init__(fig_size=fig_size)

        if subplot_kw is None:
            subplot_kw = {}

        self._axes = self._fig.add_subplot(111, **subplot_kw)

    @property
    def axes(self):
        return self._axes


class FigureWithSeveralAxes(_Figure):
    def __init__(self, n_cols, n_rows, fig_size=None):
        super().__init__(fig_size=fig_size)
        axess = self._fig.subplots(nrows=n_rows, ncols=n_cols)

        if n_rows == 1 and n_cols == 1:
            axess = [[axess]]
        elif n_rows == 1:
            axess = [list(axess)]
        elif n_cols == 1:
            axess = [[axes] for axes in axess]
        else:
            axess = [list(axess_in_row) for axess_in_row in axess]

        self._axess = axess

    def get_axes(self, row_idx, col_idx):
        return self._axess[row_idx][col_idx]


class SharedXFigure(_Figure):
    def __init__(self, num_axes, fig_size=None, sharey=False, height_ratios=None):
        super().__init__(fig_size=fig_size)
        axess = self._fig.subplots(
            nrows=num_axes,
            ncols=1,
            sharex=True,
            sharey=sharey,
            height_ratios=height_ratios,
        )
        try:
            axess[0]
        except TypeError:
            axess = [axess]
        self._axess = axess

    def get_axes(self, idx):
        return self._axess[idx]


class GeneralFigure(_Figure):
    def __init__(self, fig_size=None):
        super().__init__(fig_size=fig_size)
        self.axess = {}

    def add_axes(
        self,
        axes_id=None,
        row_idx=0,
        col_idx=0,
        left_margin=0.1,
        right_margin=0.02,
        top_margin=0.02,
        bottom_margin=0.07,
        axes_col_widths=None,
        axes_row_heights=None,
        projection=None,
        zorder=1,
        share_y=None,
    ):
        if axes_col_widths is None:
            axes_col_widths = [1]

        if axes_row_heights is None:
            axes_row_heights = [1]

        if isinstance(row_idx, int):
            row_idx = slice(row_idx, row_idx + 1)
        if isinstance(col_idx, int):
            col_idx = slice(col_idx, col_idx + 1)

        axes_col_widths = [width / sum(axes_col_widths) for width in axes_col_widths]
        axes_row_heights = [
            height / sum(axes_row_heights) for height in axes_row_heights
        ]
        assert math.isclose(sum(axes_col_widths), 1)
        assert math.isclose(sum(axes_row_heights), 1)

        this_axes_width = sum(axes_col_widths[col_idx])
        prev_cols_width = sum(axes_col_widths[: col_idx.start])
        next_cols_width = sum(axes_col_widths[col_idx.stop :])

        this_axes_height = sum(axes_row_heights[row_idx])
        prev_rows_height = sum(axes_row_heights[: row_idx.start])
        next_rows_height = sum(axes_row_heights[row_idx.stop :])

        left_margin = left_margin * this_axes_width
        right_margin = right_margin * this_axes_width
        top_margin = top_margin * this_axes_height
        bottom_margin = bottom_margin * this_axes_height

        left = prev_cols_width + left_margin
        bottom = next_rows_height + bottom_margin
        width = this_axes_width - left_margin - right_margin
        height = this_axes_height - top_margin - bottom_margin

        axes = self._fig.add_axes(
            (left, bottom, width, height),
            projection=projection,
            zorder=zorder,
            sharey=share_y,
        )

        if axes_id is None:
            axes_id = len(self.axess)
        self.axess[axes_id] = axes

        return axes


def plot_hist(values, axes, n_bins=20, range_=None):
    if range_ is None:
        range_ = numpy.nanmin(values), numpy.nanmax(values)

    bin_edges = numpy.linspace(range_[0], range_[1], num=n_bins)

    axes.hist(values, bins=bin_edges)


def plot_hist_fig(values, plot_path, n_bins=20, range_=None):
    fig = SimpleFigure()

    plot_hist(values, fig.axes, n_bins=n_bins, range_=range_)

    fig.save_fig(plot_path)


def set_y_label(axes, label):
    axes.set_ylabel(label)


def set_x_label(axes, label):
    axes.set_xlabel(label)


def set_x_ticks(axes, x_poss, labels, rotation=45, ha="right"):
    axes.set_xticklabels(labels, rotation=rotation, ha=ha)
    axes.set_xticks(x_poss)


def set_x_ticks_format(axes, rotation=45, ha="right"):
    for tick in axes.get_xticklabels():
        tick.set_rotation(rotation)
        tick.set_ha(ha)


def set_y_axis_log(axes):
    axes.set_yscale("log")


def plot_2d_joinplot(
    xvalues,
    yvalues,
    hex_cmap="viridis",
    hex_map_log=False,
    x_lim=None,
    y_lim=None,
    x_hist_n_bins=20,
    y_hist_n_bins=20,
):
    fig = GeneralFigure()

    hist_axes_width = 0.2
    hex_axes_width = 1 - hist_axes_width

    axes_col_widths = [hex_axes_width, hist_axes_width]
    axes_row_heights = [hist_axes_width, hex_axes_width]
    hex_axes = fig.add_axes(
        axes_id="hex",
        row_idx=1,
        col_idx=0,
        axes_col_widths=axes_col_widths,
        axes_row_heights=axes_row_heights,
    )
    xvalues_hist_axes = fig.add_axes(
        axes_id="xvalues_hist",
        row_idx=0,
        col_idx=0,
        axes_col_widths=axes_col_widths,
        axes_row_heights=axes_row_heights,
    )
    yvalues_hist_axes = fig.add_axes(
        axes_id="yvalues_hist",
        row_idx=1,
        col_idx=1,
        axes_col_widths=axes_col_widths,
        axes_row_heights=axes_row_heights,
    )
    bins = "log" if hex_map_log else None

    if x_lim is not None or y_lim is not None:
        mask = None
        if x_lim is not None:
            mask_x = numpy.logical_and(xvalues >= x_lim[0], xvalues <= x_lim[1])
            mask = mask_x

        if y_lim is not None:
            mask_y = numpy.logical_and(yvalues >= y_lim[0], yvalues <= y_lim[1])
            if mask is None:
                mask = mask_y
            else:
                mask = numpy.logical_and(mask, mask_y)
        xvalues = xvalues[mask]
        yvalues = yvalues[mask]

    hex_axes.hexbin(xvalues, yvalues, gridsize=20, cmap=hex_cmap, bins=bins)
    xvalues_hist_axes.hist(xvalues, bins=x_hist_n_bins)
    yvalues_hist_axes.hist(yvalues, orientation="horizontal", bins=y_hist_n_bins)
    return fig


def draw_rectangle(axes, x, y, width, height):
    rectangle = matplotlib.patches.Rectangle((x, y), width, height)
    axes.add_patch(rectangle)


def turn_off_axis(axis):
    axis.set_major_locator(matplotlib.ticker.NullLocator())


def plot_values_along_genome(
    values,
    chroms,
    poss,
    coord_converter,
    axes=None,
    plot_path=None,
    marker=".",
    linestyle="None",
    color="blue",
    chroms_are_sorted=False,
    use_scatter=False,
    marker_size=10,
    draw_chrom_lines=True,
    zorder=2,
):
    if axes is None and plot_path is None:
        raise ValueError("An axes or a plot path should be provided")
    if axes is not None and plot_path is not None:
        raise ValueError("Either axes or a plot path should be provided")

    if axes is None:
        fig = Figure((10, 5))
        FigureCanvas(fig)  # Don't remove it or savefig will fail later
        axes = fig.add_subplot(111)

    poss = numpy.array(poss)
    chroms = numpy.array(chroms)

    if not chroms_are_sorted:
        index_to_sort = chroms.argsort()
        chroms = chroms[index_to_sort]
        poss = poss[index_to_sort]
        values = values[index_to_sort]

    unique_chroms = numpy.unique(chroms)
    converted_poss = None

    for chrom in unique_chroms:
        mask = chroms == chrom
        poss_for_this_chrom = poss[mask]
        poss_for_this_chrom = coord_converter.transform_coordinate(
            chrom, poss_for_this_chrom
        )

        if converted_poss is None:
            converted_poss = poss_for_this_chrom
        else:
            converted_poss = numpy.concatenate((converted_poss, poss_for_this_chrom))

    if use_scatter:
        axes.scatter(
            converted_poss,
            values,
            marker=marker,
            color=color,
            s=marker_size,
            zorder=zorder,
        )
    else:
        axes.plot(
            converted_poss,
            values,
            marker=marker,
            linestyle=linestyle,
            color=color,
            markersize=marker_size,
            zorder=zorder,
        )

    if draw_chrom_lines:
        _draw_chromosome_lines(axes, coord_converter)

    if plot_path:
        fig.tight_layout()
        fig.savefig(str(plot_path))


def _draw_chromosome_lines(axes, coord_converter, min_y=None, max_y=None):
    y_lims = axes.get_ylim()

    if min_y is None:
        min_y = y_lims[0]
    if max_y is None:
        max_y = y_lims[1]

    try:
        chrom_lens = coord_converter.chrom_lens
    except AttributeError:
        chrom_lens = None

    if chrom_lens:
        for chrom, length in chrom_lens.items():
            chrom_end = coord_converter.transform_coordinate(chrom, 1) + length
            axes.plot(
                [chrom_end, chrom_end],
                [min_y, max_y],
                c=DARK_GRAY,
                linestyle="solid",
            )

    try:
        pericentromeric_starts = coord_converter.pericentromeric_starts
    except AttributeError:
        pericentromeric_starts = None

    if pericentromeric_starts:
        for chrom, position in pericentromeric_starts.items():
            transformed_position = coord_converter.transform_coordinate(chrom, position)
            axes.plot(
                [transformed_position, transformed_position],
                [min_y, max_y],
                c=DARK_GRAY,
            )

    xticks = []
    xtick_labels = []
    for chrom, span in coord_converter.chrom_spans.items():
        xticks.append(span[1])
        xtick_labels.append(chrom)
    axes.set_xticks(xticks)
    axes.set_xticklabels(xtick_labels, rotation=45, ha="right")

    axes.grid(False)


def _set_y_ticks(tick_poss, tick_labels, axes, rotation=0, va="center", fontsize=10):
    axes.set_yticks(tick_poss)
    axes.set_yticklabels(tick_labels, rotation=rotation, va=va, fontsize=fontsize)


def _set_x_ticks(tick_poss, tick_labels, axes, rotation=0, ha="right", fontsize=10):
    axes.set_xticks(tick_poss)
    axes.set_xticklabels(tick_labels, rotation=rotation, ha=ha, fontsize=fontsize)


def set_x_ticks(tick_poss, tick_labels, axes, rotation=0, ha="right", fontsize=10):
    _set_x_ticks(
        tick_poss, tick_labels, axes, rotation=rotation, ha=ha, fontsize=fontsize
    )


def _set_axes_background(
    axes, color="white", spine_left_color="grey", spine_bottom_color="grey"
):
    axes.set_facecolor(color)
    axes.spines["left"].set_color(spine_left_color)
    axes.spines["bottom"].set_color(spine_bottom_color)


def plot_table_classification_comparison(
    data: pandas.DataFrame,
    x_col_name: str,
    y_col_name: str,
    axes,
    color="black",
    size_multiplier=1,
    classes1_order=None,
    classes2_order=None,
):
    classification_series1 = data[x_col_name]
    classification_series2 = data[y_col_name]

    common_items = sorted(
        set(classification_series1.index).intersection(classification_series2.index),
        key=str,
    )
    classification_series1 = classification_series1.reindex(common_items)
    classification_series2 = classification_series2.reindex(common_items)

    counts = Counter(zip(classification_series1.values, classification_series2.values))

    if classes1_order is None:
        key = str
    else:
        key = lambda x: str(classes1_order.index(x)) if x in classes1_order else str(x)
    classes1 = sorted(set(key[0] for key in counts.keys()), key=key)

    if classes2_order is None:
        key = str
    else:
        key = lambda x: str(classes2_order.index(x)) if x in classes2_order else str(x)
    classes2 = sorted(set(key[1] for key in counts.keys()), key=key)

    x_poss = numpy.arange(len(classes1))
    y_poss = numpy.arange(len(classes2))

    for x_pos, class1 in zip(x_poss, classes1):
        for y_pos, class2 in zip(y_poss, classes2):
            size = counts.get((class1, class2), None)
            if size is None:
                continue
            axes.scatter(
                [x_pos], [y_pos], s=size * size_multiplier, color=color, marker="s"
            )
    _set_x_ticks(x_poss, classes1, axes, rotation=45)
    _set_y_ticks(y_poss, classes2, axes)

    _set_axes_background(axes)

    axes.set_xlabel(x_col_name)
    axes.set_ylabel(y_col_name)


def create_legend(color_for_labels, axes):
    patches = []
    for label, color in color_for_labels.items():
        patch = mpatches.Patch(color=color, label=label)
        patches.append(patch)
    axes.legend(handles=patches)


def plot_dframe_as_stacked_bars(dframe, axes, style=None, x_poss=None):
    if style is None:
        style = {}
    bar_labels = list(dframe.index)
    num_bars = len(bar_labels)
    if x_poss is None:
        x_poss = numpy.arange(num_bars)
    bottoms = numpy.zeros(num_bars)

    patches = []
    labels = []
    for bar_slice_id, values in dframe.items():
        color = style.get(bar_slice_id, {}).get("color")
        bar_container = axes.bar(
            x_poss,
            bottom=bottoms,
            height=values,
            label=bar_slice_id,
            color=color,
        )
        bottoms += values
        patches.append(bar_container.patches[0])
        labels.append(bar_slice_id)
    x_poss = pandas.Series(x_poss, index=bar_labels)

    return {"x_poss": x_poss, "patches": patches, "labels": labels}


def plot_scatter(
    xs: numpy.ndarray,
    ys: numpy.ndarray,
    axes,
    groups: numpy.ndarray = None,
    styles_for_groups: dict = None,
    plot_order: list = None,
):
    xs = numpy.array(xs)
    ys = numpy.array(ys)

    if groups is None:
        groups = numpy.zeros((xs.shape[0],), dtype=int)
        grouping = False
    else:
        groups = numpy.array(groups)
        grouping = True

    if styles_for_groups is None:
        styles_for_groups = {}

    different_groups = numpy.unique(groups)

    if plot_order:
        plot_order = {group: idx for idx, group in enumerate(plot_order)}
        different_groups = sorted(different_groups, key=lambda x: plot_order[x])

    for group in different_groups:
        style = styles_for_groups.get(group, {})
        color = style.get("color", None)
        alpha = style.get("alpha", 1)
        size = style.get("size", matplotlib.rcParams["lines.markersize"] ** 2)

        mask = groups == group
        this_xs = xs[mask]
        this_ys = ys[mask]
        label = group if grouping else None
        axes.scatter(this_xs, this_ys, label=label, color=color, alpha=alpha, s=size)


def set_spine_color(axes, color):
    for spine in axes.spines.values():
        spine.set_edgecolor(color)


def set_background_color(axes, color):
    axes.set_facecolor(color)
