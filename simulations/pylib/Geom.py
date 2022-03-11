import pylib.mix as mix
import pylib.Global_variables as GLO
import numpy as np
import matplotlib.lines as mlines


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(GLO)


class Geom:
    name = ''
    geom_type = 'NONE'
    color = 'b'
    style = '-'
    width = GLO.LINE_WIDTH * GLO.DEF_COEF_WIDTH_GEOM

    def draw(self, mpl, ax, oo):
        return


class HLine(Geom):
    geom_type = 'HLINE'
    ys = None

    def draw(self, mpl, ax, oo):
        nys = len(self.ys)
        for i in range(nys):
            mpl.axhline(
                self.ys[i],
                color=self.color,
                linestyle=self.style,
                linewidth=self.width
            )


class Line(Geom):
    geom_type = 'LINE'
    lines = []

    def draw(self, mpl, ax, oo):
        n_lines = len(self.lines)

        xs = []
        ys = []
        for i in range(n_lines):
            one_line = self.lines[i]
            xs.append([one_line[0][0], one_line[1][0]])
            ys.append([one_line[0][1], one_line[1][1]])

        for i in range(n_lines):
            one_line = mlines.Line2D(
                xs[i], ys[i],
                color=self.color,
                linestyle=self.style,
                linewidth=self.width
            )
            one_line.set_dashes([0.6, 0.6])
            ax.add_line(one_line)


class Curve(Geom):
    geom_type = 'CURVE'
    xs = None
    ys = None

    def add_curve(self, xs1, ys1):
        # xs1 - 1d array of x-coordinates
        # ys1 - 1d array of y-coordinates
        if self.xs is None:
            self.xs = []
            self.ys = []
        self.xs.append(xs1)
        self.ys.append(ys1)

    def draw(self, mpl, ax, oo):
        n_curves = len(self.ys)

        for i in range(n_curves):
            one_line = mlines.Line2D(
                self.xs[i], self.ys[i],
                color=self.color,
                linestyle=self.style,
                linewidth=self.width
            )
            one_line.set_dashes([0.6, 0.6])
            ax.add_line(one_line)


class Fill(Geom):
    geom_type = 'FILL'

    xs = None  # x-coords of nodes of a polygon (counterclockwise or clockwise)
    ys = None  # y-coords of nodes of a polygon (counterclockwise or clockwise)
    alpha = 0.2

    def draw(self, mpl, ax, oo):
        xlims = ax.get_xlim()
        ylims = ax.get_ylim()

        xs_plot = np.zeros(np.size(self.xs))
        ys_plot = np.zeros(np.size(self.xs))
        for ix in range(np.size(self.xs)):
            if self.xs[ix] == 'liml':
                xs_plot[ix] = xlims[0]
            elif self.xs[ix] == 'limr':
                xs_plot[ix] = xlims[-1]
            else:
                xs_plot[ix] = self.xs[ix]

            if self.ys[ix] == 'limb':
                ys_plot[ix] = ylims[0]
            elif self.ys[ix] == 'limu':
                ys_plot[ix] = ylims[-1]
            else:
                ys_plot[ix] = self.ys[ix]

        ax.fill(xs_plot, ys_plot, self.color, alpha=self.alpha)
        ax.set_xlim(xlims)
        ax.set_ylim(ylims)


class Annotate(Geom):
    point_to = []
    point_text = []
    line = ''
    arrowstyle = '->'
    linestyle = ':'

    def __init__(self, oo):
        self.init_from_oo(oo)

    def init_from_oo(self, oo):
        self.point_to = oo.get('point_to', [None, None])
        self.point_text = oo.get('point_text', [None, None])
        self.line = oo.get('line', '')
        self.arrowstyle = oo.get('arrowstyle', '->')
        self.linestyle = oo.get('linestyle', ':')
        self.color = oo.get('color', 'black')
        self.width = oo.get('width', 2)

    def draw(self, mpl, ax, oo):
        mpl.annotate(
            s=self.line,
            xy=(self.point_to[0], self.point_to[-1]),
            xytext=(self.point_text[0], self.point_text[-1]),
            arrowprops=dict(
                 color=self.color,
                 arrowstyle=self.arrowstyle,
                 linestyle=self.linestyle,
                 linewidth=self.width,
             )
        )


# Passing-trapped boundary in velocity domains:
def pass_trap_boundary(dd, mu):
    iar = dd['a0'] / dd['R0']
    upper_branch = np.sqrt(2 * iar * mu)
    lower_branch = - np.sqrt(2 * iar * mu)

    cone_pt = np.concatenate(
        (np.flipud(lower_branch), upper_branch)
    )
    mu_cone = np.concatenate(
        (np.flipud(mu), mu)
    )

    obj_geom = Curve()
    obj_geom.add_curve(mu_cone, cone_pt)
    obj_geom.color = 'white'

    res = {
        'geom': obj_geom,
        'lower_branch': lower_branch,
        'upper_branch': upper_branch
    }
    return res


# Parabola x = a * y^2 + x0, where x - xaxis of a plot, y - yaxis of a plot
def parabola_x(x, x0, x1, y1):
    # so far, only for positive a

    # find a:
    a = (x1 - x0) / y1 ** 2

    # work x area
    ids_x_pb, x_pb, _ = mix.get_ids(x, [x0, x[-1]])

    # new x area
    x_pb = np.concatenate((
        np.full(ids_x_pb[0], np.nan), x_pb
    ))

    # find two branches of parabola:
    upper_branch =   np.sqrt(1. / a * (x_pb - x0))
    lower_branch = - np.sqrt(1. / a * (x_pb - x0))

    cone_pb = np.concatenate(
        (np.flipud(lower_branch), upper_branch)
    )
    x_cone_pb = np.concatenate(
        (np.flipud(x_pb), x_pb)
    )

    return cone_pb, x_cone_pb, x_pb, lower_branch, upper_branch


# Vertical parallelogram, which has straight vertical sides
def parallelogram_v(x, point_left, point_right, vlength, flag_upper_points):
    # define which points are given
    sign_dir = +1
    point_ref = point_right
    if flag_upper_points:
        sign_dir = -1
        point_ref = point_left

    # upper right point:
    point_ol_y = point_left[1]  + sign_dir * vlength
    point_or_y = point_right[1] + sign_dir * vlength

    # side length:
    L = np.sqrt(
        (point_right[0] - point_left[0])**2 +
        (point_right[1] - point_left[1])**2
    )
    hlength = point_right[0] - point_left[0]
    cosTheta = np.sqrt(1 - (hlength / L) ** 2)

    # work x domain
    ids_x_pa, x_pa, _ = mix.get_ids(x, [point_left[0], point_right[0]])

    # new x area
    x_new = np.full(len(x), np.nan)
    x_new[ids_x_pa] = x_pa

    # find horizontal sides:
    dx = np.abs(x_new - point_ref[0])
    L1 = L * dx / hlength
    upper_side_y = point_ref[1] - L1 * cosTheta
    lower_side_y = upper_side_y + sign_dir * vlength

    if not flag_upper_points:
        temp_side = upper_side_y
        upper_side_y = lower_side_y
        lower_side_y = temp_side

    # point of corners:
    corners = [
        point_left,                  point_right,
        [point_left[0], point_ol_y], [point_right[0], point_or_y],
    ]

    return lower_side_y, upper_side_y, corners


# Create lines:
def lines(x, lines):
    # does not work for vertical lines

    # lines = [line, line, ...]
    # line = [point_left, point_right]
    # point_ = [x, y]

    lines_y = []
    for line_one in lines:
        point_left = line_one[0]
        point_right = line_one[1]

        # work x axis
        ids_x_line, x_line, _ = mix.get_ids(x, [point_left[0], point_right[0]])
        x_new = np.full(len(x), np.nan)
        x_new[ids_x_line] = x_line

        # side length:
        L = np.sqrt(
            (point_right[0] - point_left[0]) ** 2 +
            (point_right[1] - point_left[1]) ** 2
        )
        hlength = point_right[0] - point_left[0]
        cosTheta = np.sqrt(1 - (hlength / L) ** 2)

        # define reference point
        point_ref = point_left
        if point_left[1] < point_right[1]:
            point_ref = point_right

        # find horizontal sides:
        dx = np.abs(x_new - point_ref[0])
        L1 = L * dx / hlength
        line_y = point_ref[1] - L1 * cosTheta

        # results
        lines_y.append(line_y)
    return lines_y


# Create boundaries of an area:
def create_boundaries(nx, separate_boundaries):
    nb = len(separate_boundaries)
    boundaries = np.zeros([nx, nb])
    for ix in range(nx):
        for ib in range(nb):
            boundaries[ix, ib] = separate_boundaries[ib][ix]
    return boundaries


# Knowing boundaries of the (x,y) area, find signal in this area
def build_area(x, y, signal, boundaries):
    # boundaries: list of len(x), every element of a list is an array or list
    # of y points, which indicate area boundaries for a given x

    # - length of the axis x -
    nx = len(x)

    # - find indices for the area -
    ids_area = []
    for ix in range(nx):
        x1_boundaries = boundaries[ix]
        id_left_boundary = 0
        ids_whole_zone = np.array([])
        count_boundary = 0
        for i_boundary in range(len(x1_boundaries)):
            if np.isnan(x1_boundaries[i_boundary]):
                continue
            id_one_boundary, _, _ = mix.get_ids(y, x1_boundaries[i_boundary])
            if np.mod(count_boundary, 2) == 0:
                count_boundary += 1
                id_left_boundary = id_one_boundary
            else:
                count_boundary += 1
                ids_part_zone = np.array([i for i in range(id_left_boundary, id_one_boundary + 1)])
                ids_whole_zone = np.concatenate((
                    ids_whole_zone, ids_part_zone
                ))
        ids_area.append(ids_whole_zone)

    # - build J*E inside the strap -
    signal_area = np.full([len(x), len(y)], np.nan)
    for ix in range(nx):
        ids_x1_area = np.array(ids_area[ix]).astype(int)
        signal_area[ix, ids_x1_area] = signal[ix, ids_x1_area]

    return signal_area, ids_area



