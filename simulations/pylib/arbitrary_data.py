import pylib.mix as mix
import pylib.ControlPlot as cpr
import pylib.ymath as ymath
import numpy as np


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(cpr)
    mix.reload_module(ymath)


def choose_one_var_xy(one_signal):
    data = None if one_signal['data'] is None \
        else np.array(one_signal['data'])
    res = {
        'data': data,  # 2d data
        'x': np.array(one_signal['x']),
        'y': np.array(one_signal['y']),
        'tit': one_signal.get('title', ''),
        'x1_name': one_signal.get('x1_name', 'x'),
        'x2_name': one_signal.get('x2_name', 'y'),
        'x1_format': one_signal.get('x1_format', '{:0.3e}'),
        'x2_format': one_signal.get('x2_format', '{:0.3e}'),
        'x1_label': one_signal.get('x1_label', 'x'),
        'x2_label': one_signal.get('x2_label', 'y'),
    }
    return res


def choose_one_var_x(one_signal):
    data = None if one_signal['data'] is None \
        else np.array(one_signal['data'])

    if 'x_err' in one_signal:
        x_err = None if one_signal['x_err'] is None \
            else np.array(one_signal['x_err'])
    else:
        x_err = None

    if 'y_err' in one_signal:
        y_err = None if one_signal['y_err'] is None \
            else np.array(one_signal['y_err'])
    else:
        y_err = None

    res = {
        'data': data,  # 1d data
        'x': np.array(one_signal['x']),
        'x_err': x_err,
        'y_err': y_err,
        'tit': one_signal.get('title', ''),
        'x1_name': one_signal.get('x1_name', 'x'),
        'x1_format': one_signal.get('x1_format', '{:0.3e}'),
        'x1_label': one_signal.get('x1_label', 'x'),
    }
    return res