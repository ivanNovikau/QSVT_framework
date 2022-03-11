import pylib.mix as mix
import pylib.ymath as ymath
import pylib.curve as crv
import pylib.Global_variables as GLO
import matplotlib.pyplot as mpl
import numpy as np
import types
from matplotlib import animation, ticker
from IPython.display import HTML
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors
from matplotlib.figure import Figure


def reload():
    mix.reload_module(mix)
    mix.reload_module(crv)
    mix.reload_module(ymath)
    mix.reload_module(GLO)


def plot_curves_mix(curves):
    fig, axs = mpl.subplots(
        ncols=curves.ncols, nrows=curves.nrows,
        figsize=(curves.ff['figure_width'],
            curves.ff['figure_height']
        )
    )
    for id_col, list_curves in enumerate(curves.lists_sub_curves):
        for id_row, sub_curves in enumerate(list_curves):
            if len(list_curves) == 1:  # single row
                ax_res = axs[id_col]
            elif len(curves.lists_sub_curves) == 1:  # single column
                ax_res = axs[id_row]
            else:
                ax_res = axs[id_col, id_row]

            if sub_curves.list_curves[0].zs is None:
                plot_curves(sub_curves, fig, ax_res)
            else:
                plot_curves_3d(sub_curves, fig, ax_res)


def plot(curves, fig=None, ax=None, FIG_W=None, FIG_H=None):
    flag_1d = False
    if curves.list_curves[0].zs is None:
        flag_1d = True

    if flag_1d:
        fig, ax, css = plot_curves(curves, fig, ax, FIG_W, FIG_H)
    else:
        fig, ax, css = plot_curves_3d(curves, fig, ax, FIG_W, FIG_H)

    return fig, ax, css


def plot_curves(curves, fig=None, ax=None, FIG_W=None, FIG_H=None):
    if curves.is_empty() and not curves.flag_subplots:
        return

    figure_width = curves.ff['figure_width'] \
        if FIG_W is None else FIG_W
    figure_height = curves.ff['figure_height'] \
        if FIG_H is None else FIG_H

    # Build plots
    if not curves.ff['flag_ivis']:
        if fig is None:
            fig, ax = mpl.subplots(
                ncols=curves.ncols,
                nrows=curves.nrows,
                figsize=(figure_width, figure_height),
            )
    else:
        if fig is None:
            fig = Figure(
                figsize=(figure_width, figure_height),
            )
            ax = fig.add_subplot(111)
        else:
            ax.cla()

    # set curves
    if curves.flag_subplots:
        for id_col, list_curves in enumerate(curves.lists_sub_curves):
            for id_row, sub_curves in enumerate(list_curves):
                if len(list_curves) == 1:  # single row
                    ax_res = ax[id_col]
                elif len(curves.lists_sub_curves) == 1:  # single column
                    ax_res = ax[id_row]
                else:
                    ax_res = ax[id_col, id_row]
                set_curves(sub_curves, ax_res)
                format_plot(fig, ax_res, sub_curves)
    else:
        set_curves(curves, ax)
        format_plot(fig, ax, curves)

    return fig, ax, None


def plot_curves_3d(curves, fig=None, axs=None, FIG_W=None, FIG_H=None):
    def set_colorbar(curves, ref=None, axs=None, cax=None):
        if curves.ff['flag_colorbar']:
            boundaries_map = None
            if ref is None:
                res_colormap = GLO.DEF_COLORMAP
                first_curve = curves.list_curves[0]
                if first_curve.ff['colormap'] is not None:
                    res_colormap = first_curve.ff['colormap']

                m = mpl.cm.ScalarMappable(cmap=res_colormap)
                m.set_array(first_curve.zs.T)
                m.set_clim(curves.ff['vmin'], curves.ff['vmax'])
                ref = m
                boundaries_map = np.linspace(curves.ff['vmin'], curves.ff['vmax'], 31)

            cb = fig.colorbar(
                ref, shrink=0.8, 
                extend='both',
                cax=cax, ax=axs,
                boundaries=boundaries_map
            )

            cb.formatter.set_scientific(True)
            cb.formatter.set_powerlimits((0, 0))
            cb.ax.tick_params(
                labelsize=curves.ff['fontS'] * GLO.SCALE_TICKS
            )
            cb.ax.yaxis.get_offset_text().set_fontsize(
                curves.ff['fontS'] * GLO.SCALE_ORDER
            )

            register_offset(cb.ax.yaxis, bottom_offset)
            cb.update_ticks()

            if curves.ff['flag_graphic']:
                cb.remove()

    if curves.is_empty() and not curves.flag_subplots:
        return

    N_COLUMNS = curves.ncols
    N_ROWS = curves.nrows

    figure_width = curves.ff['figure_width'] \
        if FIG_W is None else FIG_W
    figure_height = curves.ff['figure_height'] \
        if FIG_H is None else FIG_H

    # Build plots
    if not curves.ff['flag_ivis']:
        if fig is None:
            fig, axs = mpl.subplots(
                ncols=N_COLUMNS,
                nrows=N_ROWS,
                figsize=(figure_width, figure_height),
            )
    else:
        if fig is None:
            fig = Figure(
                figsize=(figure_width, figure_height),
            )
            axs = fig.add_subplot(111)
        else:
            axs.cla()

    # set curves
    if curves.flag_subplots:
        css_res = [None for _ in range(curves.nrows * curves.ncols)]

        axs_res = [None for _ in range(curves.nrows * curves.ncols)]
        curves_res = [None for _ in range(curves.nrows * curves.ncols)]

        count_subplot = -1

        for id_col, list_curves in enumerate(curves.lists_sub_curves):
            for id_row, sub_curves in enumerate(list_curves):
                if len(list_curves) == 1:  # single row
                    ax_res = axs[id_col]
                elif len(curves.lists_sub_curves) == 1:  # single column
                    ax_res = axs[id_row]
                else:
                    ax_res = axs[id_row, id_col]

                count_subplot = count_subplot + 1

                axs_res[count_subplot]    = ax_res
                curves_res[count_subplot] = sub_curves

                css_res[count_subplot] = plot_curves_3d_subplot(sub_curves, ax_res, fig)
    else:
        css_res = plot_curves_3d_subplot(curves, axs, fig)

        css_res = [css_res]
        axs_res = [axs]
        curves_res = [curves]

    # set a colorbar for every subplot
    if curves.sel_colorbar_subplots == 'none':
        for id_cs, cs_one in enumerate(css_res):
            divider = make_axes_locatable(axs_res[id_cs])
            cax = divider.append_axes("right", "5%", pad="3%")
            set_colorbar(curves_res[id_cs], ref=cs_one, cax=cax)

    elif curves.sel_colorbar_subplots == 'all':
        divider = make_axes_locatable(mpl.gca())
        cax = divider.append_axes("right", "5%", pad="3%")
        set_colorbar(curves_res[curves.id_ref_subplot], cax=cax)

    elif curves.sel_colorbar_subplots == 'row':
        axs_row = [None] * N_COLUMNS
        for id_row in range(N_ROWS):
            for id_col in range(N_COLUMNS):
                axs_row[id_col] = axs_res[id_row + N_ROWS*id_col]

            divider = make_axes_locatable(axs_res[id_row + N_ROWS*(N_COLUMNS-1)])
            cax = divider.append_axes("right", "5%", pad="3%")

            set_colorbar(
                curves_res[id_row + N_ROWS * curves.id_ref_subplot],
                cax=cax
            )

    else:
        mix.error_mes('Wrong selector for subplot colorbar arrangement')

    # format the plot
    for id_ax, ax in enumerate(axs_res):
        set_curves(curves_res[id_ax], ax, 1)
        format_plot(fig, ax, curves_res[id_ax], flag_2d=True)

    return fig, axs, css_res


def plot_curves_3d_subplot(curves, ax, fig):
    # data from the first curve, that has to be 3d plot
    curve_one = curves.list_curves[0]
    ZZ = curve_one.zs
    if curve_one.xs.ndim < 2:
        XX, YY = np.meshgrid(curve_one.xs, curve_one.ys)
    else:
        XX = curve_one.xs
        YY = curve_one.ys

    # check colormap:
    res_colormap = GLO.DEF_COLORMAP
    if curve_one.ff['colormap'] is not None:
        res_colormap = curve_one.ff['colormap']

    # --- contour plot ---
    divnorm = None
    if curve_one.ff['colormap_center'] is not None:
        divnorm = matplotlib.colors.DivergingNorm(vcenter=curve_one.ff['colormap_center'])
    cs = ax.contourf(
        XX, YY, ZZ.T,
         levels=curve_one.ff['levels'],  # Nlevels or np.linspace(vmin, vmax, Nlevels)
         cmap=res_colormap,
         norm=divnorm,
         vmin=curves.ff['vmin'],  # for proper change of vmax, change levels as well as np.linspace(vmin, vmax, Nlevels)
         vmax=curves.ff['vmax'],
         extend='both'
    )
    return cs


def set_curves(curves, ax, id_curve_start=0):
    for icrv in range(id_curve_start, curves.n()):
        curve = curves.list_curves[icrv]

        y_res = curve.ys
        if y_res is None:
            continue

        # legend: list to one line:
        res_legend = mix.create_line_from_list(curve.ff['legend'])

        # check color:
        res_color = curve.ff['color'] \
            if curve.ff['color'] is not None \
            else GLO.new_color(icrv)

        # check style:
        res_style = curve.ff['style']
        if res_style is None:
            res_style = GLO.new_style(icrv) \
                if curves.ff['flag_diff_styles'] \
                else GLO.DEF_ONE_STYLE

        # create a curve on a plot
        if curve.ff['flag_hist']:
            if res_legend is not None:
                ax.hist(y_res, curve.xs, alpha=curve.ff['pr_alpha'],
                        label=res_legend, density=True)
            else:
                ax.hist(y_res, curve.xs, alpha=curve.ff['pr_alpha'],
                        density=True)
        else:
            if curves.ff['flag_norm']:
                y_res = ymath.find_norm(y_res, curve.ff['norm_to'])
            if curves.ff['flag_semilogy']:
                y_res = abs(y_res)

            if not curve.ff['flag_errorbar']:
                if curves.ff['flag_semilogy']:
                    ref_lines, = ax.semilogy(curve.xs, abs(y_res), res_style)
                else:
                    ref_lines, = ax.plot(curve.xs, y_res, res_style)
                ref_line_format = ref_lines
            else:
                ref_lines = ax.errorbar(curve.xs, curve.ys,
                                        yerr=curve.ys_err,
                                        xerr=curve.xs_err,
                                        fmt=res_style,
                                        elinewidth=curve.ff['width'] * GLO.ERR_WIDTH_COEF,
                                        ecolor=res_color)
                ref_line_format = ref_lines[0]

            # set line and marker sizes
            if res_style == ':':
                ref_line_format.set_dashes(curve.ff['dashed_format'])
            mpl.setp(ref_line_format,
                     linewidth=curve.ff['width'],
                     color=res_color,
                     markersize=curve.ff['markersize'],
                     markerfacecolor=curve.ff['markerfacecolor'],
                     markeredgewidth=curve.ff['width'] * GLO.MARKER_EDGE_WIDTH_COEF)

            # set legend
            if res_legend is not None:
                ref_lines.set_label(res_legend)


def format_plot(fig, ax, curves, flag_2d=False):
    ncurves = curves.n()
    ngeoms = curves.n_geoms
    ntexts = len(curves.list_text)

    # set fixed limits
    if curves.ff['flag_fixed_limits']:
        curves.set_fixed_limits()

    # set labels:
    res_xlabel = mix.create_line_from_list(curves.ff['xlabel'])
    res_ylabel = mix.create_line_from_list(curves.ff['ylabel'])

    if res_xlabel is None:
        res_xlabel = ""
    if res_ylabel is None:
        res_ylabel = ""

    ax.set_xlabel(
        res_xlabel,
        fontsize=curves.ff['fontS'] * GLO.SCALE_LABELS
    )
    ax.set_ylabel(
        res_ylabel,
        fontsize=curves.ff['fontS'] * GLO.SCALE_LABELS
    )

    # axes ticks:
    if not curves.ff['flag_ivis']:
        mpl.sca(ax)

    sci_limits = curves.ff['sci_limits']

    if not np.isnan(curves.ff['xticks']).any():
        ax.set_xticks(curves.ff['xticks'])
    if not np.isnan(curves.ff['xticks_labels']).any():
        ax.set_xticklabels(curves.ff['xticks_labels'])
    else:
        ax.ticklabel_format(axis='x', style=curves.ff['x_style'], scilimits=sci_limits)
        if curves.ff['flag_maxlocator']:
            ax.xaxis.set_major_locator(mpl.MaxNLocator(curves.ff['maxlocator']))

    if not np.isnan(curves.ff['yticks']).any():
        ax.set_yticks(curves.ff['yticks'])
    if not np.isnan(curves.ff['yticks_labels']).any():
        ax.set_yticklabels(curves.ff['yticks_labels'])
    else:
        if not curves.ff['flag_semilogy']:
            ax.ticklabel_format(axis='y', style=curves.ff['y_style'], scilimits=sci_limits)

    # fontsize of axes ticks
    ax.xaxis.set_tick_params(
        labelsize=curves.ff['fontS'] * GLO.SCALE_TICKS
    )
    ax.yaxis.set_tick_params(
        labelsize=curves.ff['fontS'] * GLO.SCALE_TICKS
    )
    ax.xaxis.get_offset_text().set_fontsize(
        curves.ff['fontS'] * GLO.SCALE_ORDER
    )
    ax.yaxis.get_offset_text().set_fontsize(
        curves.ff['fontS'] * GLO.SCALE_ORDER
    )
    register_offset(ax.yaxis, top_offset)

    # set limits:
    if curves.ff['xlimits'] is not None:
        ax.set_xlim(curves.ff['xlimits'][0], curves.ff['xlimits'][-1])
    if curves.ff['ylimits'] is not None:
        ax.set_ylim(curves.ff['ylimits'][0], curves.ff['ylimits'][-1])

    # set legend
    flag_legend_res = curves.ff['flag_legend']
    legends = []
    for one_curve in curves.list_curves:
        legends.append(one_curve.ff['legend'])
    if all(leg is None for leg in legends):
        flag_legend_res = False

    if ncurves > 1 and flag_legend_res:
        ax.legend(
            fontsize=curves.ff['fontS'] * GLO.LEG_SCALE,
            loc=curves.ff['legend_position'],
            facecolor=curves.ff['legend_fcol'],
            labelspacing=0.1,
            handlelength=1,
            handletextpad=0.4
        )

    # set title
    res_title = mix.create_line_from_list(curves.ff['title'])

    if res_title is None or curves.ff['flag_graphic']:
        res_title = ""

    ax.set_title(res_title,
              fontsize=curves.ff['fontS'] * GLO.SCALE_TITLE,
              pad=curves.ff['pad_title'],
              usetex=True)
    if GLO.FLAG_LATEX:
        mpl.rc('text', usetex=True)
        # mpl.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
        # mpl.rcParams['text.latex.preamble'] = [r'\boldmath']
        mpl.rcParams['text.latex.preamble'] = r"\usepackage{amsmath} \boldmath"

    # draw geometrical figures:
    for igeom in range(ngeoms):
        one_geom = curves.list_geoms[igeom]
        if one_geom is not None:
            one_geom.draw(mpl, ax, {})

    # add text:
    ax.texts = []
    if curves.ff['flag_add_text_plot']:
        for itext in range(ntexts):
            loc_text = curves.list_text[itext]
            ax.text(
                loc_text.x,
                loc_text.y,
                loc_text.line,
                fontsize=GLO.FONT_SIZE * loc_text.coef_width,
                color=loc_text.color,
                ha="center",
            )
            ax.texts[-1].set_visible(not loc_text.flag_invisible)

    # set grid
    if not flag_2d:
        ax.grid(True)

    # keep only graphic:
    if curves.ff['flag_graphic']:
        frame1 = fig.gca()
        frame1.axes.get_xaxis().set_visible(False)
        frame1.axes.get_yaxis().set_visible(False)

    if GLO.FLAG_LATEX and curves.ff['flag_tight_layout']:
        fig.tight_layout()


def bottom_offset(self, bboxes, bboxes2):
    bottom = self.axes.bbox.ymin
    self.offsetText.set(va="top", ha="left")
    self.offsetText.set_position(
        (0, bottom - self.OFFSETTEXTPAD * 8 * self.figure.dpi / 72.0))


def top_offset(self, bboxes, bboxes2):
    top = self.axes.bbox.ymax
    self.offsetText.set(va="bottom", ha="left")
    self.offsetText.set_position(
        (0, top + self.OFFSETTEXTPAD * 6 * self.figure.dpi / 72.0))


def register_offset(axis, func):
    axis._update_offset_text_position = types.MethodType(func, axis)


def animation_1d(curves):
    # WORKS, but in development

    # PLOTTING along Y, ANIMATION along X, DATA is Z

    # number of curves
    ncurves = curves.n()
    ngeoms = curves.n_geoms

    # Build plots
    fig, ax = mpl.subplots(figsize=(GLO.FIG_SIZE_W, GLO.FIG_SIZE_H))
    axes = mpl.gca()

    # set limits:
    if curves.xlimits is not None:
        ax.set_xlim(curves.xlimits[0], curves.xlimits[-1])
    if curves.ylimits is not None:
        ax.set_ylim(curves.ylimits[0], curves.ylimits[-1])

    # set labels:
    if curves.xlabel is not None:
        mpl.xlabel(r'\boldmath $' + curves.xlabel + '$', fontsize=GLO.FONT_SIZE * 1.7)
    if curves.ylabel is not None:
        mpl.ylabel(r'\boldmath $' + curves.ylabel + '$', fontsize=GLO.FONT_SIZE * 1.7)

    # axes ticks:
    if curves.xticks_labels is np.nan:
        mpl.xticks(curves.xticks) if curves.xticks is not np.nan else 0
    else:
        mpl.xticks(curves.xticks, curves.xticks_labels) if curves.xticks is not np.nan else 0

    if curves.yticks_labels is np.nan:
        mpl.yticks(curves.yticks) if curves.yticks is not np.nan else 0
    else:
        mpl.yticks(curves.yticks, curves.yticks_labels) if curves.yticks is not np.nan else 0

    # fontsize of axes ticks
    ax.xaxis.set_tick_params(labelsize=GLO.FONT_SIZE)
    ax.yaxis.set_tick_params(labelsize=GLO.FONT_SIZE)
    ax.xaxis.get_offset_text().set_fontsize(GLO.FONT_SIZE)
    ax.yaxis.get_offset_text().set_fontsize(GLO.FONT_SIZE)

    # format of axis labels
    mpl.ticklabel_format(axis='x', style=curves.x_style, scilimits=(-2, 2))
    if curves.flag_maxlocator:
        ax.xaxis.set_major_locator(mpl.MaxNLocator(curves.maxlocator))

    if curves.flag_semilogy is False:
        mpl.ticklabel_format(axis='y', style=curves.y_style, scilimits=(-2, 2))

    # set legend
    if curves.flag_legend:
        ax.legend(fontsize=GLO.FONT_SIZE * GLO.LEG_SCALE, loc=curves.legend_position,
                  facecolor=curves.legend_fcol)

    # set title
    if curves.title is not None:
        mpl.title(r'\boldmath $' + curves.title + '$', fontsize=GLO.FONT_SIZE * 1.5, pad='18', usetex=True)
    if GLO.FLAG_LATEX:
        mpl.rc('text', usetex=True)
        mpl.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
        mpl.rcParams['text.latex.preamble'] = [r'\boldmath']

    # set grid
    mpl.grid(True)

    # set empty curves
    ref_lines = [None] * ncurves
    for icrv in range(ncurves):
        curve = curves.list(icrv)
        if not curve.flag_errorbar:
            if curves.flag_semilogy:
                ref_lines[icrv], = ax.semilogy([], [], curve.style)
            else:
                ref_lines[icrv], = ax.plot([], [], curve.style)

            if curve.style == ':':
                ref_lines[icrv].set_dashes([0.5, 0.4])
            mpl.setp(ref_lines[icrv], linewidth=curve.width,
                     color=curve.color,
                     markersize=curve.markersize,
                     markerfacecolor=curve.markerfacecolor,
                     markeredgewidth=curve.width / 2)
        else:
            ref_lines[icrv] = ax.errorbar([], [],
                                    yerr=curve.ys_err, xerr=curve.xs_err, fmt=curve.style,
                                    elinewidth=curve.width / 3, ecolor=curve.color)
            mpl.setp(ref_lines[icrv][0], linewidth=curve.width,
                     color=curve.color,
                     markersize=curve.markersize,
                     markerfacecolor=curve.markerfacecolor,
                     markeredgewidth=curve.width / 2)

        # set legend
        if curve.legend == "_":
            ref_lines[icrv].set_label("_")
        else:
            ref_lines[icrv].set_label(r'\boldmath $' + curve.legend + '$')

    # # draw geometrical figures:
    # for igeom in range(ngeoms):
    #     one_geom = curves.list_geoms[igeom]
    #     one_geom.draw(mpl, ax, axes, {})

    # if FLAG_LATEX:
    #     fig.tight_layout()

    # initialization function: plot the background of each frame
    def init():
        for icrvL in range(ncurves):
            ref_lines[icrvL].set_data([], [])
        return ref_lines,  # !!!

    # animation function. This is called sequentially
    def animate(i, Y_res, Z_res):
        for icrvL in range(ncurves):
            if i < np.shape(z_res)[0]:
                ref_lines[icrvL].set_data(Y_res[icrvL][:], Z_res[icrvL][i, :])
        return ref_lines,  # !!!

    nx_max = np.zeros(ncurves)
    for icrv in range(ncurves):
        nx_max[icrv] = np.size(curves.list(icrv).xs)
    nx_max = int(np.max(nx_max))

    Y_res, Z_res = [], []
    for icrv in range(ncurves):
        curve = curves.list(icrv)

        z_res = curve.zs
        if z_res is None:
            continue

        if curves.flag_norm:
            z_res = ymath.find_norm(z_res, curve.data_norm_to)
        if curves.flag_semilogy:
            z_res = abs(z_res)
        Z_res.append(z_res)
        Y_res.append(curve.ys)

    anim = animation.FuncAnimation(
        fig,
        animate, fargs=(Y_res, Z_res),
        init_func=init,
        frames=nx_max, interval=20, blit=True
    )

    HTML(anim.to_html5_video())