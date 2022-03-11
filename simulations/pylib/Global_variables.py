from IPython.core.getipython import get_ipython
import numpy as np


def reload():
    return


def to_rgb(rgb):
    # rgb = (int, int, int)
    return "#%02x%02x%02x" % rgb

dd_null = {'project_name': 'NULL\ PROJECT'}

# type compatible with the complex-like type that ORB5 understands
comp_datatype = np.dtype([
    ('real', np.float),
    ('imaginary', np.float)
])

# ---------------------------------------------------------------------------
MIN_N_PEAKS = 3
COEF_ERR = 1.96          # alpha-quantile of the standard normal distribution
CONFIDENCE_PERC = 0.95   # (2*alpha - 1)-confidence interval, that corresponds to the alpha quantile

# --- System float ---
if isinstance(0.1, np.float64):
    SYS_FLOAT_TYPE = np.float64
else:
    SYS_FLOAT_TYPE = np.float32

# --- DEFAULT POST-PROCESSING OPERATIONS ---
NONE_FILTER = {'sel_filt': None}
DEF_FILTER_SMOOTH = {
    'operation': 'filtering',
    'domain': None,
    'oo_filters': [
        {'sel_filt': 'smooth', 'norm_w': 1, 'wind': 3}
    ]
}
DEF_FILTER_ROUGH = {
    'operation': 'filtering',
    'domain': None,
    'oo_filters': [
        {'sel_filt': 'rough', 'norm_w': 1, 'w_interval': [0, 1e-3]}
    ]
}
DEF_FFT = {
    'flag_f2': False
}
DEF_FFT_2D = {
    'flag_f2': False,
    'name_coord_fft': 't',
}
DEF_OPERATION_FFT_1D = {
    'operation': 'fft-1d',
    'oo_fft': DEF_FFT,
}
DEF_OPERATION_FFT_2D = {
    'operation': 'fft-2d',
    'oo_fft': DEF_FFT_2D,
}

# --- FOR PLOTTING ---
if 'Terminal' in get_ipython().__class__.__name__:
    # FLAG_LATEX = True
    FLAG_LATEX = False
    # FLAG_IVIS = True
    FLAG_IVIS = False
else:
    FLAG_LATEX = False
    FLAG_IVIS = False

DASHES_FORMAT = [0.5, 0.4]
MARKER_EDGE_WIDTH_COEF = 0.5
ERR_WIDTH_COEF = 0.33
DEF_MAXLOCATOR = 6
COLORMAP_LEVELS = 60
DEF_ONE_COLOR = 'blue'
DEF_ONE_STYLE = '-'
DEF_COLORMAP = 'jet'
DEF_COEF_WIDTH_GEOM = 0.5
DEF_NORM_WG = ''
DEF_COLORS = ['b', 'r', 'g', 'black', 'm', 'c',
              'y', 'k', 'cyan', 'Purple', 'gray', 'lightcoral']
DEF_STYLES = ['-', ':', '-.', ':']
DEF_SIGN_M = 1
DEF_TITLE_PAD = 18
if FLAG_LATEX:
    if FLAG_IVIS:
        FLAG_LATEX = True
        FIG_SIZE_W = 15
        FIG_SIZE_H = 9.5
        LEG_SCALE = 1.5
        FONT_SIZE = 28
        SCALE_LABELS = 1.8
        SCALE_TICKS  = 1.5
        SCALE_ORDER  = 1.2
        SCALE_TITLE  = 1.3
        LINE_WIDTH = 6
        MARKER_SIZE = 14
    else:
        FLAG_LATEX = True
        FIG_SIZE_W = 30
        FIG_SIZE_H = 19
        LEG_SCALE = 1
        FONT_SIZE = 22
        SCALE_LABELS = 1.8
        SCALE_TICKS = 1.5
        SCALE_ORDER = 1.2
        SCALE_TITLE = 1.3
        LINE_WIDTH = 6
        MARKER_SIZE = 14
else:
    FLAG_LATEX = False
    # FIG_SIZE_W = 10
    # FIG_SIZE_H = 6
    FIG_SIZE_W = 16
    FIG_SIZE_H = 10
    LEG_SCALE = 0.5
    FONT_SIZE = 22
    SCALE_LABELS = 0.6
    SCALE_TICKS = 0.6
    SCALE_ORDER = 0.6
    SCALE_TITLE = 0.6
    LINE_WIDTH = 6
    MARKER_SIZE = 14
DEF_PLOT_FORMAT = {  # describe format of a plot
    'xlabel': None,
    'ylabel': None,
    'zlabel': None,
    'wlabel': None,
    'title': None,
    'flag_semilogy': False,
    'flag_norm': False,
    'flag_colorbar': True,
    'fontS': FONT_SIZE,
    'xlimits': None,
    'ylimits': None,
    'zlimits': None,
    'xticks': np.nan,
    'yticks': np.nan,
    'ivis_add_xticks': np.nan,
    'ivis_add_yticks': np.nan,
    'xticks_labels': np.nan,
    'yticks_labels': np.nan,
    'flag_legend': True,
    'legend_position': 'best',  # 'upper right', 'center left', 'lower left'
    'legend_fcol': 'lightgray',
    'flag_diff_styles': False,
    'x_style': 'sci',  # 'sci', 'plain'
    'y_style': 'sci',  # 'sci', 'plain'
    'flag_maxlocator': False,
    'maxlocator': DEF_MAXLOCATOR,
    'flag_fixed_limits': False,
    'figure_width': FIG_SIZE_W,
    'figure_height': FIG_SIZE_H,
    'pad_title': DEF_TITLE_PAD,
    'vmin': None,
    'vmax': None,
    'flag_tight_layout': True,
    'sci_limits': (-2, 2),
    'flag_graphic': False,
    'flag_ivis': FLAG_IVIS,
    'flag_add_text_plot': True,  # to render or not an additional text
}
DEF_CURVE_FORMAT = {  # describe format of a curve
    'legend': None,
    'style': None,
    'width': LINE_WIDTH,
    'color': None,
    'markersize': MARKER_SIZE,
    'markerfacecolor': "None",
    'colormap': DEF_COLORMAP,  # hot, jet, pink, hot_r, jet_r etc.
    'colormap_center': None,
    'levels': COLORMAP_LEVELS,
    'pr_alpha': 1,
    'flag_errorbar': False,
    'flag_hist': False,
    'norm_to': None,
    'dashed_format': DASHES_FORMAT,
}


def new_color(count_color):
    if count_color < len(DEF_COLORS):
        one_color = DEF_COLORS[count_color]
    else:
        one_color = ','.join('{}'.format(*k) for k in enumerate(np.random.rand(3, 1)))
        one_color = 'rgb({})'.format(one_color)
    return one_color


def new_style(count_style):
    if count_style < len(DEF_STYLES):
        one_style = DEF_STYLES[count_style]
    else:
        one_style = ':'
    return one_style

# --- oracletool ---
OT_label_color = to_rgb((200, 200, 200))
OT_frame_color = to_rgb((120, 120, 120))
OT_selected_element_inf_frame = to_rgb((168, 172, 185))
OT_canvas_color = to_rgb((140, 140, 140))
OT_border_color = to_rgb((60, 60, 60))
OT_color_tabs_frame = to_rgb((160, 160, 160))
OT_color_button = to_rgb((160, 160, 160))
OT_color_active_button = to_rgb((203, 203, 255))
OT_COEF_FONT_SIZE = 3
OT_height_tab_frame = 30


# --- DEFAULT VARIABLE DEFINITIONS ---
DEF_SPECIES = 'deuterium'
def_erbar_ts = {
    'type':             'zonal',
    'variable':         'er',
    'plane':            'ts',
    'avr_operation':    'point-s',
    'avr_domain':       0.5,
}
def_potsc_chi1_t = {
    'type': 'nonzonal',
    'variable': 'potsc',
    'plane': 'ts',
    'chi-point': 0.0,
    'avr_operation': 'point-s',
    'avr_domain': 0.5,
}
def_fields3d_ts = {
    'type': 'fields3d',
    'plane': 'ts',
}
def_fields3d_n1 = {
    'type': 'fields3d',
    'variable': 'n1',
    'plane': 'ts',
    'n1': np.nan,
    'chi-point': 0.0,
}
def_fields3d_n1_schi = {
    'type': 'fields3d',
    'variable': 'n1',
    'plane': 'schi',
    'n1': np.nan,
    'chi-point': 0.0,
}
def_potsc_rz = {
    'type': 'nonzonal',
    'variable': 'potsc',
    'plane': 'rz',
    't-point': 0.0,
    'phi-point': 0.0,
}
def_safety_factor = {
    'type': 'equ-profile',
    'variable': 'q',
    'plane': 'ts',
    'avr_operation': 'point-t',
    'avr_domain': 0,
}
def_Teq_deuterium = {
    'type': 'equ-profile',
    'variable': 'T-keV',
    'species': 'deuterium',
    'plane': 'ts',
    'avr_operation': 'point-t',
    'avr_domain': 0,
}
def_n_deuterium = {
    'type': 'equ-profile',
    'variable': 'n',
    'species': 'deuterium',
    'plane': 'ts',
    'avr_operation': 'point-t',
    'avr_domain': 0,
}
def_T_evol_deuterium = {
    'type': 'transport',
    'variable': 'T',
    'species': 'deuterium',
    'plane': 'ts',
}
def_T_keV_evol_deuterium = {
    'type': 'transport',
    'variable': 'T-keV',
    'species': 'deuterium',
    'plane': 'ts',
}
def_je = {
    'type':             'mpr',
    'variable':         'je',
    'species':          DEF_SPECIES,
    'mu-domain':        None,
    'vpar-domain':      None,
    'flag-VparMuIntegration-ComplexDomain': False,
    'ids-vpar-area': None,
    'line-name-area': '',
}
def_je_t = {
    'type':             'mpr',
    'variable':         'je',
    'plane':            'tnone',
    'avr_operation':    'none-',
    'species':          DEF_SPECIES,
    'mu-domain':        None,
    'vpar-domain':      None,
    'flag-VparMuIntegration-ComplexDomain': False,
    'ids-vpar-area':    None,
    'line-name-area': '',
}
def_efield = {
    'type':             'mpr',
    'variable':         'efield',
    'plane':            'tnone',
    'avr_operation':    'none-',
    'species':          'total',
}
def_chi0_ts = {
    'type':             'transport',
    'variable':         'chi_norm0',
    'plane':            'ts',
    'avr_operation':    'point-s',
    'avr_domain':       0.5,
}
def_arbitrary_1d = {
    'type':             'arbitrary',
    'plane':            'xnone',
    'avr_operation':    'none-',
}
def_arbitrary_2d = {
    'type':             'arbitrary',
    'plane':            'xy',
}


def create_signal(default_signal, dd):
    res_signal = dict(default_signal)
    res_signal.update({
        'dd': dd
    })
    return res_signal


def create_signals_dds(default_signal, dds,
                      types=None, variables=None, species=None,
                      planes=None, operations=None, domains=None,
                      flag_arbitrary=False, xs=None, ys=None, datas=None,
                       xs_err=None, ys_err=None):
    n_signals = len(dds)
    res_signals = []
    for id_signal in range(n_signals):
        one_type        = get_field(id_signal, types,       default_signal.get('type', None))
        one_variable    = get_field(id_signal, variables,   default_signal.get('variable', None))
        one_species     = get_field(id_signal, species,     DEF_SPECIES)
        one_plane       = get_field(id_signal, planes,      default_signal.get('plane', None))
        one_operation   = get_field(id_signal, operations,  default_signal.get('avr_operation', None))
        one_domain      = get_field(id_signal, domains,     default_signal.get('avr_domain', None))

        one_signal = dict(default_signal)
        one_signal.update({
            'dd': dds[id_signal],
            'type': one_type,
            'variable': one_variable,
            'species': one_species,
            'plane': one_plane,
            'avr_operation': one_operation,
            'avr_domain': one_domain,
        })

        if flag_arbitrary:
            x    = get_field(id_signal, xs, None)
            y    = get_field(id_signal, ys, None)
            x_err = get_field(id_signal, xs_err, None)
            y_err = get_field(id_signal, ys_err, None)
            data = get_field(id_signal, datas, None)
            one_signal.update({
                'x': x,
                'y': y,
                'data': data,
                'x_err': x_err,
                'y_err': y_err,
            })

        res_signals.append(one_signal)
    return res_signals


def get_field(id_field, fields, default_field):
    if fields is None:
        return default_field
    return fields[id_field] \
        if id_field < len(fields) \
        else default_field


# set values_field to fields with names names_field
#  to a dictionary or array of dictionaries struc
def set_fields(struc, names_field, values_field):
    # for every dictionary in list struc_list
    # set one value from values_field
    # to one field from names_field
    if isinstance(struc, list):
        struc_list = list(struc)
    else:
        struc_list = [struc]

    for id_dict, one_dict in enumerate(struc_list):
        one_dict[names_field[id_dict]] = values_field[id_dict]


