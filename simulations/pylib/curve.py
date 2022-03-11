import pylib.mix as mix
import pylib.Global_variables as GLO

import numpy as np

def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(GLO)

class Curve:
    CurveName = ''
    ff = None  # format of the curve
    xs = None
    xs_err = None
    ys_err = None
    ys = None
    zs = None
    ws = None

    def __init__(self):
        self.ff = dict(GLO.DEF_CURVE_FORMAT)

    def name(self, v):
        self.CurveName = v
        return self

    def XS(self, v):
        self.xs = v
        return self

    def YS(self, v):
        self.ys = v
        return self

    def ZS(self, v):
        self.zs = v
        return self

    def WS(self, v):
        self.ws = v
        return self

    def set_ff(self, v):
        self.ff = dict(v)
        return self

    def set_errorbar(self, v, ys=None, xs=None):
        self.ff['flag_errorbar'] = v
        self.ys_err = ys
        self.xs_err = xs
        return self

    def copy(self):
        res = Curve()
        res.xs = np.array(self.xs) if self.xs is not None else None
        res.ys = np.array(self.ys) if self.ys is not None else None
        res.zs = np.array(self.zs) if self.zs is not None else None

        res.xs_err = np.array(self.xs_err) if self.xs_err is not None else None
        res.ys_err = np.array(self.ys_err) if self.ys_err is not None else None
        res.ws = np.array(self.ws) if self.ws is not None else None

        res.ff = dict(self.ff)
        res.CurveName = self.CurveName

        return res

    def get_flag_2d(self):
        if self.zs is None:
            return False
        else:
            return True


class Curves:
    list_curves = None
    map_curves = None
    n_curves = 0
    list_text = None
    list_geoms = None
    n_geoms = 0
    ff = None  # format

    lists_sub_curves = None     # to plot curves in different subplots:
                                # every element - one object of class Curves
                                # 2D list: [id_col, id_row]
    ncols, id_col_set = 1, 0
    nrows, id_row_set = 1, 0
    flag_subplots = False
    sel_colorbar_subplots = 'none'  # 'row', 'column', 'all'
    id_ref_subplot = 0  # if 'all' - subplot with id_ref_subplot will define a common colobar
                     # if 'row' - in every row, subplot with id_ref_subplot will define a common colobar
                     # if 'column' - in every column, subplot with id_ref_subplot will define a common colobar

    def __init__(self):
        self.list_curves = []
        self.map_curves = {}

        self.list_geoms = []
        self.list_text = []

        # create default format:
        self.ff = GLO.DEF_PLOT_FORMAT

    def new(self, name_curve=None):
        new_curve = Curve()
        self.list_curves.append(new_curve)
        self.n_curves += 1
        if name_curve is None:
            name_curve = 'curve_' + str(self.n_curves-1)
        self.map_curves[name_curve] = new_curve.name(name_curve)

        return new_curve

    def create_sub(self, ncols=1, nrows=1,
                   selector_colorbar_subplots='none', id_reference_subplot=0):
        self.flag_subplots = True
        self.ncols = ncols
        self.nrows = nrows
        self.lists_sub_curves = [
            [None for _ in range(self.nrows)]
            for _ in range(self.ncols)
        ]
        self.sel_colorbar_subplots = selector_colorbar_subplots
        self.id_ref_subplot = id_reference_subplot

    def put_sub(self, curves_to_put, id_col=None, id_row=None):
        id_col_res = id_col if id_col is not None else self.id_col_set
        id_row_res = id_row if id_row is not None else self.id_row_set
        self.lists_sub_curves[id_col_res][id_row_res] = curves_to_put

        # define new identifiers for the sub_curves list:
        id_col_loc = -1
        for el in self.lists_sub_curves:
            id_col_loc += 1
            if None in el:
                self.id_col_set = id_col_loc
                self.id_row_set = el.index(None)
                break

    def load(self, curves_to_load):
        if curves_to_load is None:
            return
        for one_curve in curves_to_load.list_curves:
            self.n_curves += 1
            self.list_curves.append(one_curve)
            self.map_curves[one_curve.name] = one_curve

    def set_ff(self, v):
        self.ff = dict(v)
        return self

    def newg(self, geoms):
        if isinstance(geoms, list):
            for one_geom in geoms:
                self.list_geoms.append(one_geom)
                self.n_geoms += 1
        else:
            self.list_geoms.append(geoms)
            self.n_geoms += 1

    def newt(self, oo_text):
        if isinstance(oo_text, list):
            for one_oo_text in oo_text:
                self.list_text.append(PlText(one_oo_text))
        else:
            self.list_text.append(PlText(oo_text))

    def n(self):
        return self.n_curves

    def set_fixed_limits(self):
        cr1 = self.list_curves[0]
        if cr1.xs is not None:
            self.ff['xlimits'] = [cr1.xs[0], cr1.xs[-1]]
        if cr1.ys is not None:
            self.ff['ylimits'] = [cr1.ys[0], cr1.ys[-1]]
        if cr1.zs is not None:
            self.ff['zlimits'] = [cr1.zs[0], cr1.zs[-1]]

    def sort_legends(self):
        count_curve = -1
        new_list_curves = []
        list_curves_emp_leg = []
        for one_curve in self.list_curves:
            count_curve = count_curve + 1
            leg = one_curve.ff['legend']

            if leg is None:
                list_curves_emp_leg.append(one_curve)
            else:
                new_list_curves.append(one_curve)
        new_list_curves.extend(list_curves_emp_leg)
        self.list_curves = new_list_curves

    def is_empty(self):
        if len(self.list_curves) == 0:
            return True
        else:
            return False

    def get_legends(self):
        legends = []
        for id_curve, ocurve in enumerate(self.list_curves):
            leg = ocurve.ff['legend']
            if leg is None:
                leg = "curve {:d}".format(id_curve)
            legends.append(leg)

        return legends


class PlText:
    x = None
    y = None
    line = ''
    color = 'black'
    coef_width = 1
    flag_invisible = False
    oo_init = None

    def __init__(self, oo):
        self.oo_init = oo
        self.init_from_oo(oo)

    def init_from_oo(self, oo):
        self.x = oo.get('x', None)  # in units of x-axis
        self.y = oo.get('y', None)  # in units of y-axis
        self.color = oo.get('color', 'black')
        self.coef_width = oo.get('coef_width', 1)
        self.line = mix.create_line_from_list(oo.get('line', ''))


def copy_curves(curves, ax):
    curves_copy = Curves()

    # copy format
    curves_copy.set_ff(curves.ff)

    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    xticks, yticks = ax.get_xticks(), ax.get_yticks()
    xticks = xticks[1:]  if xlim[0] > xticks[0] else xticks
    xticks = xticks[:-1] if xlim[-1] < xticks[-1] else xticks
    yticks = yticks[1:] if ylim[0] > yticks[0] else yticks
    yticks = yticks[:-1] if ylim[-1] < yticks[-1] else yticks

    # curves_copy.ff.update({
    #     'xticks': xticks,
    #     'yticks': yticks,
    # })

    # copy curves
    for one_curve in curves.list_curves:
        new_curve = one_curve.copy()

        curves_copy.n_curves += 1
        curves_copy.list_curves.append(new_curve)
        curves_copy.map_curves[one_curve.name] = new_curve

    # copy texts
    for one_text in curves.list_text:
        curves_copy.list_text.append(
            PlText(one_text.oo_init)
        )

    # copy geometries
    curves_copy.newg(curves.list_geoms)

    return curves_copy






