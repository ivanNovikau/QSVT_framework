import pylib.mix as mix
import numpy as np
import scipy.signal
from scipy.stats import norm as stat_norm

import pylib.Global_variables as GLO
import pylib.ymath as ymath
import pylib.arbitrary_data as ARD
import pylib.Geom as geom
import pylib.curve as crv
import pylib.ControlPlot as cpr


def reload():
    mix.reload_module(mix)
    mix.reload_module(GLO)
    mix.reload_module(ymath)
    mix.reload_module(ARD)
    mix.reload_module(geom)
    mix.reload_module(crv)
    mix.reload_module(cpr)
    

def choose_vars(oo):
    def x1x2_format(vv, xx1, xx2):
        vv['x1'], vv['fx1'], vv['labx'] = xx1[0], xx1[1], xx1[2]
        vv['x2'], vv['fx2'], vv['laby'] = xx2[0], xx2[1], xx2[2]

    if 'signals' not in oo:
        oo_signals = [oo.get('signal', None)]
    else:
        oo_signals = oo.get('signals', [])

    count_signal, vvars = -1, []
    for one_signal in oo_signals:
        count_signal += 1

        # choose type of the signal
        opt_type = one_signal['type']
        ref_module = None
        # if opt_type == 'zonal':
        #     ref_module = zf
        # elif opt_type == 'transport':
        #     ref_module = transport
        # elif opt_type == 'nonzonal':
        #     ref_module = itg
        # elif opt_type == 'equ-profile':
        #     ref_module = equil_profiles
        # elif opt_type == 'fields3d':
        #     ref_module = fields3d
        # elif opt_type.lower() == 'distribution':
        #     ref_module = distribution
        # elif opt_type.lower() == 'mpr':
        #     ref_module = MPR
        if opt_type.lower() == 'arbitrary':
            ref_module = ARD
        else:
            mix.error_mes('Wrong signal type.')

        # choose coordinate system, where the signal will be considered:
        opt_plane = one_signal['plane']
        vvar_plane = None
        if opt_plane == 'ts':
            vvar_plane = ref_module.choose_one_var_ts(one_signal)
            x1x2_format(vvar_plane,
                        ['t', '{:0.3e}', 't'],
                        ['s', '{:0.3f}', 's'])
        elif opt_plane == 'tchi':
            vvar_plane = ref_module.choose_one_var_tchi(one_signal)
            x1x2_format(vvar_plane,
                        ['t', '{:0.3e}', 't'],
                        ['chi', '{:0.3f}', '\chi'])
        elif opt_plane == 'tvpar':
            vvar_plane = ref_module.choose_one_var_tvpar(one_signal)
            x1x2_format(vvar_plane,
                        ['t', '{:0.3e}', 't'],
                        ['vpar', '{:0.3f}', 'v_{\parallel}'])
        elif opt_plane == 'tnone':
            vvar_plane = ref_module.choose_one_var_t(one_signal)
            x1x2_format(vvar_plane,
                        ['t', '{:0.3e}', 't'],
                        [None, None, None])
        elif opt_plane == 'vparmu':
            vvar_plane = ref_module.choose_one_var_vparmu(one_signal)
            x1x2_format(vvar_plane,
                        ['mu', '{:0.3f}', '\mu'],
                        ['vpar', '{:0.3f}', 'v_{\parallel}'])
        elif opt_plane == 'rz':
            vvar_plane = ref_module.choose_one_var_rz(one_signal)
            x1x2_format(vvar_plane,
                        ['r', '{:0.3f}', 'R'],
                        ['z', '{:0.3f}', 'Z'])
        elif opt_plane == 'schi':
            vvar_plane = ref_module.choose_one_var_rz(one_signal)
            x1x2_format(vvar_plane,
                        ['s', '{:0.3f}', 's'],
                        ['chi', '{:0.3f}', '\chi'])
        elif opt_plane == 'xy':
            vvar_plane = ref_module.choose_one_var_xy(one_signal)
            x1x2_format(vvar_plane,
                [vvar_plane['x1_name'], vvar_plane['x1_format'], vvar_plane['x1_label']],
                [vvar_plane['x2_name'], vvar_plane['x2_format'], vvar_plane['x2_label']],
            )
        elif opt_plane == 'xnone':
            vvar_plane = ref_module.choose_one_var_x(one_signal)
            x1x2_format(vvar_plane,
                [vvar_plane['x1_name'], vvar_plane['x1_format'], vvar_plane['x1_label']],
                [None, None, None],
            )
        else:
            mix.error_mes('Wrong name of plane.')

        # set reference signal:
        one_signal.update({'flag_var_first': True}) if count_signal == 0 else \
            one_signal.update({'flag_var_first': False})
        if one_signal['flag_var_first'] and 'x2' in vvar_plane:
            oo.update({vvar_plane['x1'] + '_ref': vvar_plane[vvar_plane['x1']]})
            if vvar_plane['x2'] is not None:
                oo.update({vvar_plane['x2'] + '_ref': vvar_plane[vvar_plane['x2']]})

        # averaging of the chosen signal
        vvar = ymath.avr_x1x2(vvar_plane, one_signal, oo)
 
        # *** signal legend ***
        # pr_name = one_signal['dd']['project_name'] if 'dd' in one_signal else ''
        # pr_name += ':\ ' if pr_name is not '' else ''
        # vvar['leg'] = pr_name + vvar['line_avr'] + ':\ ' + vvar['tit']
        vvar['leg'] = vvar['line_avr'] + ':\ ' + vvar['tit']

        # save signal
        vvars.append(vvar)
    return vvars

def build_plot(curves, fig=None, ax=None):
    # --- create a figure and plot data ---
    fig, ax, css = cpr.plot(
        curves, fig, ax,
        FIG_W=curves.ff['figure_width']/2,
        FIG_H=curves.ff['figure_height']/2,
    )
    return fig, ax, css

def plot_vars_1d(oo):
    # signals to plot
    vvars = choose_vars(oo)
    signals = oo.get('signals', [])
    n_vars = len(vvars)

    # - additional data -
    ff = dict(oo.get('ff', GLO.DEF_PLOT_FORMAT))  # format
    oo_text             = oo.get('text', [])
    geoms               = oo.get('geoms', [])
    sel_norm_x          = oo.get('sel_norm_x', None)
    sel_norm_ys         = oo.get('sel_norm_ys', [None])
    oo_postprocessing   = oo.get('oo_postprocessing', None)
    flag_plot           = oo.get('flag_plot', True)

    # normalization (first stage):
    line_x_norm = mix.normalization(sel_norm_x)['line_norm']
    line_y_norm = mix.normalization(sel_norm_ys[0])['line_norm'] \
        if len(sel_norm_ys) == 1 else ''
    if len(sel_norm_ys) == 1:
        sel_norm_ys = [sel_norm_ys[0]] * n_vars

    # XY labels
    if ff['xlabel'] is not None:
        ff['xlabel'] += line_x_norm
    if ff['ylabel'] is not None:
        ff['ylabel'] += line_y_norm

    # Create a plot
    curves = crv.Curves().set_ff(ff)

    # additional text and geometrical figures:
    curves.newt(oo_text)
    curves.newg(geoms)

    # styles, colors, legends
    stys    = ff.get('styles', [])
    colors  = ff.get('colors', [])
    legends = ff.get('legends', [])
    flags_hist = ff.get('flags_hist', None)

    # - different variables -
    for ivar in range(n_vars):
        vvar = vvars[ivar]
        data = vvar['data']
        if data is None:
            continue
        x     = np.array(vvar['x'])
        x_err = vvar.get('x_err', None)
        y_err = vvar.get('y_err', None)
        leg  = vvar['leg']
        dd_one = signals[ivar]['dd'] if 'dd' in signals[ivar] else None
        oo_var_operations = oo_postprocessing[ivar] \
            if oo_postprocessing is not None else None

        # curve format
        ff_curve = dict(GLO.DEF_CURVE_FORMAT)

        # different flags:
        if flags_hist is not None:
            ff_curve['flag_hist'] = flags_hist[ivar]

        # normalization (second stage):
        coef_x_norm = 1
        if not ff_curve['flag_hist']:
            coef_x_norm = mix.normalization(sel_norm_x, dd_one)['coef_norm']

        sel_norm_y = sel_norm_ys[ivar] if ivar < len(sel_norm_ys) else None
        temp_dict = mix.normalization(sel_norm_y, dd_one)
        line_leg_norm = temp_dict['line_norm']
        coef_y_norm   = temp_dict['coef_norm']

        # - post-processing -
        if not ff_curve['flag_hist']:
            data, x = ymath.post_processing(data, x, oo_var_operations)

        # domain of plotting:
        # add x_end, x_start in your option
        # to change limits of plots with rescaling of the plot
        if not ff_curve['flag_hist']:
            x, ids_x = mix.get_array_oo(oo, x, 'x')
            data = mix.get_slice(data, ids_x)

        # x normalization
        if not ff_curve['flag_hist']:
            x = x * coef_x_norm
        data = data * coef_y_norm

        # style, color
        ff_curve['style'] = stys[ivar]   if ivar < len(stys)   else None
        ff_curve['color'] = colors[ivar] if ivar < len(colors) else None

        # legend
        one_leg = \
            leg + [line_leg_norm] if isinstance(leg, list) else \
                leg + line_leg_norm
        ff_curve['legend'] = legends[ivar] if len(legends) > ivar else one_leg

        # - add a new curve -
        curves.new().XS(x).YS(data).set_ff(ff_curve)
        if x_err is not None or y_err is not None:
            curves.list_curves[-1].set_errorbar(True, ys=y_err, xs=x_err)

    # - plot the curves -
    if len(curves.list_curves) is not 0 and flag_plot:
        # cpr.plot_curves(curves)
        # ivis.plot_data(curves)
        build_plot(curves)

    if not flag_plot:
        return curves
    else:
        return None

def plot_several_curves(oo):
    list_curves = oo.get('list_curves', None)
    if list_curves is None:
        return
    flag_subplots = oo.get('flag_subplots', False)
    flag_3d       = oo.get('flag_3d', False)
    flag_mix      = oo.get('flag_mix', False)

    # combine all plots
    curves_result = crv.Curves()
    if not flag_subplots:
        count_element = -1
        curves_ref = None
        for current_curves in list_curves:
            count_element += 1
            if count_element == 0:
                curves_ref = current_curves
            curves_result.load(current_curves)

        # set styling:
        ff = dict(oo.get('ff', curves_ref.ff))
        curves_result.set_ff(ff)

        # styling
        for id_curve, one_curve in enumerate(curves_result.list_curves):
            ff_curve = dict(GLO.DEF_CURVE_FORMAT)
            for field_curve in ff_curve.keys():
                temp = ff.get(field_curve + 's', [])
                if len(temp) > id_curve:
                    ff_curve[field_curve] = temp[id_curve]
                else:
                    if one_curve.ff[field_curve] != ff_curve[field_curve]:
                        ff_curve[field_curve] = one_curve.ff[field_curve]
            one_curve.set_ff(ff_curve)

        # plot curves
        if len(curves_result.list_curves) is not 0:
            if not flag_3d:
                cpr.plot_curves(curves_result)
            else:
                cpr.plot_curves_3d(curves_result)
    else:
        ncols = oo.get('ncols', 1)
        nrows = oo.get('nrows', 1)
        sel_colorbar_subplots = oo.get('sel_colorbar_subplots', 'none')
        id_ref_subplot = oo.get('id_ref_subplot', 0)

        ff_global = oo.get('ff', dict(GLO.DEF_PLOT_FORMAT))
        curves_result.set_ff(ff_global)

        curves_result.create_sub(
            ncols, nrows,
            selector_colorbar_subplots=sel_colorbar_subplots,
            id_reference_subplot=id_ref_subplot
        )

        # different Curves objects from the list_curves
        # are put consequently to the subplot matrix COLUMN by COLUMN:
        count_curves = -1
        for id_col in range(ncols):
            for id_row in range(nrows):
                count_curves += 1
                if count_curves < len(list_curves):
                    curves_result.put_sub(
                        list_curves[count_curves], id_col, id_row
                    )
                else:
                    break

        if flag_mix:
            cpr.plot_curves_mix(curves_result)
        else:
            if not flag_3d:
                cpr.plot_curves(curves_result)
            else:
                cpr.plot_curves_3d(curves_result)

def calc_wg(oo, oo_wg):
    # -------------------------------------------------------------------------------
    # -> oo_var - dictionary to choose a variable
    #   (input dict. for the function choose_vars(...))
    # -------------------------------------------------------------------------------
    # -> oo_wg - dictionary with parameters to calculate frequency and dynamic rate:
    # 't_work' - work time domain
    # 'sel_wg' - line, name of a method to calculate frequency and rate:
    #   'wg-adv', 'w-adv', 'g-adv', 'wg-est', 'w-est', 'g-est'
    # 'flag_two_stages' = True:
    #       Firstly, one calculates gamma,
    #       then create a signal = initial_signal * exp(-gamma*t) and
    #       calculate the frequency
    #   False:
    #       calculate gamma and frequency from the initial_signal
    # 'sel_norm': 'wc', 'vt', 'khz':
    #       output normalization
    # ---
    # 'flag_stat' = True:
    #       calculate errorbars
    # 'n_samples' - integer:
    #       number of time interval variations
    # 'min_n_peaks' - integer:
    #       minimum number of peaks in one time interval
    # 'threshold_w' - float:
    #       relative difference between estimated (linear fitting) value of frequency
    #       (w_est) and
    #       frequency value found from NL fitting (w_adv),
    #       if |(w_adv - w_est)/w_est| <= threshold_w, then we are taking w_adv as
    #       a result frequency, otherwise we don't take any value
    # 'threshold_g' - float:
    #       the same as threshold_w, but for the damping rate
    # ---
    # - FILTERING -
    #  -> If 'flag_two_stages' = True, there are three stages of the filtering:
    #       global, for gamma, for frequency;
    #  -> Globally filtered signal is a starting signal for the calculation of the
    #       both gamma and frequency;
    #  -> After that, globally filtered signal can be filtered separately
    #       before the calculation of the gamma and before the calc. of the frequency
    #  -> If 'flag_two_stages' = False, there is only global filtering
    # 'filt_global' - dict. or [dict., dict., ...]:
    #       global filtering
    # 'filt_gamma' - dict. or [dict., dict., ...]:
    #       additional filtering of the globally filtered signal
    #       before the calculation of the gamma
    # 'filt_freq' - dict. or [dict., dict., ...]:
    #       additional filtering of the globally filtered signal
    #       before the calculation of the frequency
    #  -> For the description of these dictionaries, see the function ymath.filtering
    # -------------------------------------------------------------------------------
    # -> oo_plot - dictionary for plotting:
    # 't_plot' - domain of plotting;
    # 'flag_norm' = True: normalized plots;
    # 'flag_semilogy' = True: Y-axis in logarithmic scale;
    # 'flag_plot_print' = True: plot results and print values on screen
    # -------------------------------------------------------------------------------

    # - None-filter -
    non_filt = GLO.NONE_FILTER
    out_res = {}

    # - FUNCTION: filtering at one stage -
    def one_stage_filtering(x, y, oo_filt_loc):
        oo_filts_res = []
        if oo_filt_loc is None or len(oo_filt_loc) is 0:
            oo_filts_res.append(non_filt)
        elif isinstance(oo_filt_loc, dict):
            oo_filts_res.append(oo_filt_loc)
        else:
            oo_filts_res = oo_filt_loc  # array of filters

        # apply one by one the filters in the array :
        dict_loc = {'x': np.array(x), 'filt': np.array(y)}
        count_filt = -1
        for one_filt in oo_filts_res:
            count_filt += 1
            dict_loc = ymath.filtering(
                dict_loc['x'], dict_loc['filt'], one_filt
            )
        return dict_loc

    # - FUNCTION: get result -
    def give_res(dict_wg_loc, name_method, name_value, coef_norm):
        value_res, err_value = None, None
        if dict_wg_loc[name_method] is not None:
            value_res = dict_wg_loc[name_method][name_value]
            err_value = dict_wg_loc[name_method][name_value + '_err']

        # normalization:
        if value_res is not None:
            value_res *= coef_norm
        if err_value is not None:
            err_value *= coef_norm

        line_value = 'None'
        if value_res is not None:
            line_value = '{:0.3e}'.format(value_res)
            if err_value is not None:
                line_value += ' +- {:0.3e}'.format(err_value)

        return value_res, line_value

    # - FUNCTION: different options of w and g measurement -
    def find_wg(t, y, sel_wg, flag_print=True):
        res_dict = {}
        if sel_wg == 'wg-adv':
            res_dict = ymath.advanced_wg(t, y, flag_print=flag_print)
        elif sel_wg == 'w-adv':
            res_dict = ymath.advanced_w(t, y, flag_print=flag_print)
        elif sel_wg == 'g-adv':
            res_dict = ymath.advanced_g(t, y, flag_print=flag_print)
        elif sel_wg == 'wg-est':
            res_dict = {
                'est': ymath.estimate_wg(t, y),
                'adv': None
            }
        elif sel_wg == 'w-est':
            res_dict = {
                'est': ymath.estimate_w(t, y),
                'adv': None
            }
        elif sel_wg == 'g-est':
            res_dict = {
                'est': ymath.estimate_g(t, y),
                'adv': None
            }
        else:
            mix.error_mes('Wrong wg-selector: check oo_wg[''sel_wg'']')

        return res_dict

    # - project structure -
    signal = oo.get('signal', None)
    ff = dict(oo.get('ff', dict(GLO.DEF_PLOT_FORMAT)))
    flag_plot_print = ff.get('flag_plot_print', True)
    dd = signal.get('dd', {})

    # seperate plots or subplots:
    flag_subplots = oo.get('flag_subplots', False)
    flag_plot_internal = not flag_subplots

    # - choose a variable -
    dict_var = choose_vars(oo)[0]
    leg_data = dict_var['leg']

    # - initial data -
    data_init = dict_var['data']
    t_init    = dict_var['x']
    dict_fft_init = ymath.filtering(t_init, data_init, non_filt)

    # --- Frequency/rate calculation PARAMETERS ---
    t_work = oo_wg.get('t_work', [])
    sel_wg = oo_wg.get('sel_wg', 'wg-adv')

    if len(t_work) == 0 or t_work is None:
        t_work = t_init
    flag_two_stages = oo_wg.get('flag_two_stages', False)
    if sel_wg == 'w-adv' or sel_wg == 'w-est' or\
            sel_wg == 'g-adv' or sel_wg == 'g-est':
        flag_two_stages = False

    oo_filt_global = oo_wg.get('filt_global', None)
    oo_filt_gamma  = oo_wg.get('filt_gamma', None)
    oo_filt_freq   = oo_wg.get('filt_freq', None)

    flag_stat   = oo_wg.get('flag_stat', False)
    n_samples   = oo_wg.get('n_samples', None)
    min_n_peaks = oo_wg.get('min_n_peaks', None)
    threshold_w = oo_wg.get('threshold_w', 0.1)
    threshold_g = oo_wg.get('threshold_g', 0.2)
    n_bins      = oo_wg.get('n_bins', 40)

    flag_print = False
    if flag_stat:
        if n_samples <= 10:
            flag_print = True

    # normalization
    coef_norm_global_w, coef_norm_global_g, line_norm_w, line_norm_g = \
        mix.choose_wg_normalization(
            oo_wg.get('sel_norm_wg', GLO.DEF_NORM_WG), dd
        )

    # --- GLOBAL FILTERING ---
    dict_global = one_stage_filtering(t_init, data_init, oo_filt_global)
    data_global = dict_global['filt']
    t_global    = dict_global['x']

    # --- WORK DOMAIN ---
    ids_work, t_work, _ = mix.get_ids(t_global, t_work)
    data_work = data_global[ids_work]
    dict_fft_work_global    = ymath.filtering(t_work, data_work, non_filt)
    ids_peaks_work, _       = scipy.signal.find_peaks(np.abs(data_work))
    t_peaks_work = t_work[ids_peaks_work]

    # --- PLOTTING: GLOBAL FILTERING ---
    if flag_plot_print:
        # signal
        nsignals = 3  # original, globally filtered, absolute peaks
        ch_signals_time_evol = GLO.create_signals_dds(
            GLO.def_arbitrary_1d, [dd] * nsignals,
            flag_arbitrary=True,
            xs=[t_init, t_global, t_peaks_work],
            datas=[data_init, data_global, data_work[ids_peaks_work]],
        )

        # styling:
        ff_time_evol = dict(ff)
        ff_time_evol.update({
            'xlabel': 't',
            'title': leg_data,
            'legends': ['original', 'globally\ filtered', 'peaks'],
            'styles': ['-', ':', 'o'],
        })

        # geometry
        area_work = geom.Fill()
        area_work.xs = [t_work[0], t_work[-1], t_work[-1], t_work[0]]
        area_work.ys = ['limb', 'limb', 'limu', 'limu']
        area_work.color = 'grey'
        area_work.alpha = 0.3

        # plotting
        oo_time_evolution = {
            'signals': ch_signals_time_evol,
            'ff': ff_time_evol,
            'geoms': [area_work],
            'flag_plot': flag_plot_internal,
        }
        curves_t = plot_vars_1d(oo_time_evolution)
        del ch_signals_time_evol, ff_time_evol, oo_time_evolution

        # - FAST FOURIER TRANSFORM -
        # signal
        nsignals = 3  # original, globally filtered, globally filtered in a working domain
        ch_signals_fft = GLO.create_signals_dds(
            GLO.def_arbitrary_1d, [dd] * nsignals,
            flag_arbitrary=True,
            xs=[
                dict_fft_init['w2'],
                dict_global['w2'],
                dict_fft_work_global['w2']
            ],
            datas=[
                dict_fft_init['fft_init_2'],
                dict_global['fft_filt_2'],
                dict_fft_work_global['fft_init_2']
            ],
        )

        # styling:
        ff_fft = dict(ff)
        ff_fft.update({
            'xlabel': '\omega',
            'title': 'FFT',
            'legends': [
                'FFT:\ initial',
                ['FFT:\ globally\ filtered:', 'whole\ time\ domain'],
                ['FFT:\ globally\ filtered:', 'work\ time\ domain']],
            'styles': ['-', ':', ':'],
            'flag_semilogy': False,
        })

        # plotting
        oo_fft = {
            'signals': ch_signals_fft,
            'ff': ff_fft,
            'flag_plot': flag_plot_internal,
        }
        curves_fft = plot_vars_1d(oo_fft)
        del ch_signals_fft, ff_fft, oo_fft

        # plot subplots if necessary:
        if flag_subplots:
            oo_sub = {
                'ncols': 1,
                'nrows': 2,
                'list_curves': [curves_t, curves_fft],
                'flag_subplots': flag_subplots,
            }
            plot_several_curves(oo_sub)

    # --- NAIVE CALCULATION ---
    if not flag_two_stages:
        dict_wg = find_wg(t_work, data_work, sel_wg=sel_wg, flag_print=flag_print)

        # - results -
        w_est, line_w_est = give_res(dict_wg, 'est', 'w', coef_norm_global_w)
        g_est, line_g_est = give_res(dict_wg, 'est', 'g', coef_norm_global_g)
        w_adv, line_w_adv = give_res(dict_wg, 'adv', 'w', coef_norm_global_w)
        g_adv, line_g_adv = give_res(dict_wg, 'adv', 'g', coef_norm_global_g)

        # --- PLOT FITTING ---
        if flag_plot_print:
            # signal
            nsignals = 3  # work signal, peaks, linear regression
            flag_adv = False
            if dict_wg['adv'] is not None:
                flag_adv = True
                nsignals = 4  # + NL FITTING

            ch_signals = GLO.create_signals_dds(
                GLO.def_arbitrary_1d, [dd] * nsignals,
                flag_arbitrary=True,
                xs=[
                    t_work,
                    t_work[dict_wg['est']['ids_peaks']],
                    dict_wg['est']['x_fit']
                ] + [dict_wg['adv']['x_fit']] if flag_adv else [],
                datas=[
                    data_work,
                    data_work[dict_wg['est']['ids_peaks']],
                    dict_wg['est']['y_fit']
                ] + [dict_wg['adv']['y_fit']] if flag_adv else [],
            )

            # styling:
            ff_current = dict(ff)
            ff_current.update({
                'xlabel': 't',
                'title': 'Freq./Gamma\ calculation',
                'legends': [
                    leg_data, 'peaks', 'LIN.\ REGRESSION'
                ] + ['NL\ FITTING']  if flag_adv else [],
                'styles': [
                    '-', 'o', ':'
                ] + [':'] if flag_adv else [],
            })

            # plotting
            oo_current = {
                'signals': ch_signals,
                'ff': ff_current,
            }
            plot_vars_1d(oo_current)

            # Print results
            line_res = '--- NAIVE CALCULATION ---\n'
            line_res += '--- ESTIMATION ---\n'
            line_res += 'w' + line_norm_w + ' = ' + line_w_est + '\n'
            line_res += 'g' + line_norm_g + ' = ' + line_g_est + '\n'
            line_res += '--- NL FITTING ---\n'
            line_res += 'w' + line_norm_w + ' = ' + line_w_adv + '\n'
            line_res += 'g' + line_norm_g + ' = ' + line_g_adv
            print(line_res)
    else:
        # - FIND GAMMA -
        dict_gamma_filt = one_stage_filtering(t_work, data_work, oo_filt_gamma)
        data_gamma = dict_gamma_filt['filt']
        t_gamma    = dict_gamma_filt['x']
        dict_gamma = find_wg(t_gamma, data_gamma, sel_wg, flag_print=flag_print)

        w_est_prel, line_w_est_prel = give_res(dict_gamma, 'est', 'w', coef_norm_global_w)
        g_est, line_g_est           = give_res(dict_gamma, 'est', 'g', coef_norm_global_g)
        w_adv_prel, line_w_adv_prel = give_res(dict_gamma, 'adv', 'w', coef_norm_global_w)
        g_adv, line_g_adv           = give_res(dict_gamma, 'adv', 'g', coef_norm_global_g)

        # - FIND FREQUENCY -
        g_mult = g_est
        if g_adv is not None:
            g_mult = g_adv
        g_mult /= coef_norm_global_g
        data_work_exp = data_work * np.exp(- g_mult * t_work)
        dict_freq_filt = one_stage_filtering(
            t_work,
            data_work_exp,
            oo_filt_freq
        )
        data_freq = dict_freq_filt['filt']
        t_freq    = dict_freq_filt['x']
        dict_freq = find_wg(t_freq, data_freq, sel_wg, flag_print=flag_print)

        w_est, line_w_est           = give_res(dict_freq, 'est', 'w', coef_norm_global_w)
        g_est_zero, line_g_est_zero = give_res(dict_freq, 'est', 'g', coef_norm_global_g)
        w_adv, line_w_adv           = give_res(dict_freq, 'adv', 'w', coef_norm_global_w)
        g_adv_zero, line_g_adv_zero = give_res(dict_freq, 'adv', 'g', coef_norm_global_g)

        # --- PLOTTING ---
        if flag_plot_print:

            # --- GAMMA: TIME EVOLUTION ---
            # signal
            nsignals = 2  # original, filtered
            ch_signals = GLO.create_signals_dds(
                GLO.def_arbitrary_1d, [dd] * nsignals,
                flag_arbitrary=True,
                xs=[t_work, t_gamma],
                datas=[data_work, data_gamma],
            )

            # styling:
            ff_current = dict(ff)
            ff_current.update({
                'xlabel': 't',
                'title': 'Gamma\ calculation:\ filtering',
                'legends': ['origianl', 'filtered'],
                'styles': ['-', ':'],
            })

            # plotting
            oo_current = {
                'signals': ch_signals,
                'ff': ff_current,
                'flag_plot': flag_plot_internal,
            }
            curves_gt = plot_vars_1d(oo_current)

            # --- GAMMA: FFT ---
            # signal
            nsignals = 2  # globally filtered, filtered
            ch_signals = GLO.create_signals_dds(
                GLO.def_arbitrary_1d, [dd] * nsignals,
                flag_arbitrary=True,
                xs=[dict_gamma_filt['w2'], dict_gamma_filt['w2']],
                datas=[dict_gamma_filt['fft_init_2'], dict_gamma_filt['fft_filt_2']],
            )

            # styling:
            ff_current = dict(ff)
            ff_current.update({
                'xlabel': '\omega',
                'title': 'Gamma\ calculation:\ FFT',
                'legends': ['FFT:\ global\ filtering', 'FFT:\ gamma\ filtering'],
                'styles': ['-', ':'],
                'flag_semilogy': False,
            })

            # plotting
            oo_current = {
                'signals': ch_signals,
                'ff': ff_current,
                'flag_plot': flag_plot_internal,
            }
            curves_gfft = plot_vars_1d(oo_current)

            # --- GAMMA: FITTING ---
            # signal
            nsignals = 3  # work signal, peaks, linear regression
            flag_adv = False
            if dict_gamma['adv'] is not None:
                flag_adv = True
                nsignals = 4  # + NL FITTING

            ch_signals = GLO.create_signals_dds(
                GLO.def_arbitrary_1d, [dd] * nsignals,
                flag_arbitrary=True,
                xs=[
                       t_gamma,
                       t_work[dict_gamma['est']['ids_peaks']],
                       dict_gamma['est']['x_fit']
                   ] + [dict_gamma['adv']['x_fit']] if flag_adv else [],
                datas=[
                          data_gamma,
                          data_work[dict_gamma['est']['ids_peaks']],
                          dict_gamma['est']['y_fit']
                      ] + [dict_gamma['adv']['y_fit']] if flag_adv else [],
            )

            # styling:
            ff_current = dict(ff)
            ff_current.update({
                'xlabel': 't',
                'title': 'Gamma\ calculation:\ fitting',
                'legends': [
                               'filtered',
                               'peaks',
                               'LIN.\ REGRESSION'
                           ] + ['NL\ FITTING'] if flag_adv else [],
                'styles': [
                              '-', 'o', ':'
                          ] + [':'] if flag_adv else [],
            })

            # plotting
            oo_current = {
                'signals': ch_signals,
                'ff': ff_current,
                'flag_plot': flag_plot_internal,
            }
            curves_gfit = plot_vars_1d(oo_current)

            # plot gamma plots
            if flag_subplots:
                oo_sub = {
                    'ncols': 1,
                    'nrows': 3,
                    'list_curves': [curves_gt, curves_gfft, curves_gfit],
                    'flag_subplots': flag_subplots,
                }
                plot_several_curves(oo_sub)

            # --- FREQUENCY: TIME EVOLUTION ---
            # signal
            nsignals = 2  # original, filtered
            ch_signals = GLO.create_signals_dds(
                GLO.def_arbitrary_1d, [dd] * nsignals,
                flag_arbitrary=True,
                xs=[t_work, t_freq],
                datas=[data_work, data_freq],
            )

            # styling:
            ff_current = dict(ff)
            ff_current.update({
                'xlabel': 't',
                'title': 'Freq.\ calculation:\ filtering',
                'legends': ['original', '*\exp(-g*t),\ filtered'],
                'styles': ['-', ':'],
            })

            # plotting
            oo_current = {
                'signals': ch_signals,
                'ff': ff_current,
                'flag_plot': flag_plot_internal,
            }
            curves_wt = plot_vars_1d(oo_current)

            # --- FREQUENCY: FFT ---
            # signal
            nsignals = 2  # globally filtered, filtered
            ch_signals = GLO.create_signals_dds(
                GLO.def_arbitrary_1d, [dd] * nsignals,
                flag_arbitrary=True,
                xs=[dict_freq_filt['w2'], dict_freq_filt['w2']],
                datas=[dict_freq_filt['fft_init_2'], dict_freq_filt['fft_filt_2']],
            )

            # styling:
            ff_current = dict(ff)
            ff_current.update({
                'xlabel': '\omega',
                'title': 'Freq.\ calculation:\ FFT',
                'legends': ['FFT:\ global\ filtering', 'FFT:\ freq.\ filtering'],
                'styles': ['-', ':'],
                'flag_semilogy': False,
            })

            # plotting
            oo_current = {
                'signals': ch_signals,
                'ff': ff_current,
                'flag_plot': flag_plot_internal,
            }
            curves_wfft = plot_vars_1d(oo_current)

            # --- FREQUENCY: FITTING ---
            # signal
            nsignals = 3  # work signal, peaks, linear regression
            flag_adv = False
            if dict_freq['adv'] is not None:
                flag_adv = True
                nsignals = 4  # + NL FITTING

            ch_signals = GLO.create_signals_dds(
                GLO.def_arbitrary_1d, [dd] * nsignals,
                flag_arbitrary=True,
                xs=[
                       t_freq,
                       t_work[dict_freq['est']['ids_peaks']],
                       dict_freq['est']['x_fit']
                   ] + [dict_freq['adv']['x_fit']] if flag_adv else [],
                datas=[
                          data_freq,
                          data_work_exp[dict_freq['est']['ids_peaks']],
                          dict_freq['est']['y_fit']
                      ] + [dict_freq['adv']['y_fit']] if flag_adv else [],
            )

            # styling:
            ff_current = dict(ff)
            ff_current.update({
                'xlabel': 't',
                'title': 'Freq.\ calc.:\ fitting',
                'legends': [
                               'filtered',
                               'peaks',
                               'LIN.\ REGRESSION'
                           ] + ['NL\ FITTING'] if flag_adv else [],
                'styles': [
                              '-', 'o', ':'
                          ] + [':'] if flag_adv else [],
                'norm_to': data_work_exp
            })

            # plotting
            oo_current = {
                'signals': ch_signals,
                'ff': ff_current,
                'flag_plot': flag_plot_internal,
            }
            curves_wfit = plot_vars_1d(oo_current)

            # plot frequency plots
            if flag_subplots:
                oo_sub = {
                    'ncols': 1,
                    'nrows': 3,
                    'list_curves': [curves_wt, curves_wfft, curves_wfit],
                    'flag_subplots': flag_subplots,
                }
                plot_several_curves(oo_sub)

            # - print results -
            line_res = '--- NAIVE CALCULATION ---\n'
            line_res += '- GAMMA: ESTIMATION -\n'
            line_res += 'prel. w' + line_norm_w + ' = ' + line_w_est_prel + '\n'
            line_res += 'g' + line_norm_g + ' = ' + line_g_est + '\n'
            line_res += '- GAMMA: NL FITTING -\n'
            line_res += 'prel. w' + line_norm_w + ' = ' + line_w_adv_prel + '\n'
            line_res += 'g' + line_norm_g + ' = ' + line_g_adv + '\n'
            line_res += '- FREQUENCY: ESTIMATION -\n'
            line_res += 'w' + line_norm_w + ' = ' + line_w_est + '\n'
            line_res += '(g_real - g_num)' + line_norm_g + ' = ' + line_g_est_zero + '\n'
            line_res += '- FREQUENCY: NL FITTING -\n'
            line_res += 'w' + line_norm_w + ' = ' + line_w_adv + '\n'
            line_res += '(g_real - g_num)' + line_norm_g + ' = ' + line_g_adv_zero
            print(line_res)

    # - save results of naive calculation -
    naive_res = {
        'w_est': w_est, 'g_est': g_est,
        'w_adv': w_adv, 'g_adv': g_adv,
    }
    out_res.update({'naive': naive_res})

    # --- CALCULATION WITH STATISTICS ---
    t_intervals = []
    if flag_stat:
        oo_get_intervals = {
            'nsamples': n_samples,
            't_work': t_work,
            'min_n_periods': min_n_peaks,
            't_period': np.mean(np.diff(t_peaks_work))
        }
        dict_intervals = mix.get_t_intervals(oo_get_intervals, flag_plot_print)
        t_intervals = dict_intervals['t_intervals']

        # plot time intervals
        if flag_print and flag_plot_print:
            for one_t_interval in t_intervals:
                # signal
                nsignals = 2  # globally filtered, peaks
                ch_signals = GLO.create_signals_dds(
                    GLO.def_arbitrary_1d, [dd] * nsignals,
                    flag_arbitrary=True,
                    xs=[t_global, t_work[ids_peaks_work]],
                    datas=[data_global, data_work[ids_peaks_work]],
                )

                # styling:
                ff_current = dict(ff)
                ff_current.update({
                    'xlabel': 't[\omega_{ci}^{-1}]',
                    'title': 'Chosen\ time\ intervals',
                    'legends': [
                        'Globally\ filtered\ data', 'peaks'],
                    'styles': ['-', 'o'],
                })

                # geometry
                area_calc_chosen = geom.Fill()
                area_calc_chosen.xs = [
                    one_t_interval[0], one_t_interval[-1],
                    one_t_interval[-1], one_t_interval[0]
                ]
                area_calc_chosen.ys = ['limb', 'limb', 'limu', 'limu']
                area_calc_chosen.color = 'grey'
                area_calc_chosen.alpha = 0.6

                # plotting
                oo_current = {
                    'signals': ch_signals,
                    'ff': ff_current,
                    'geoms': [area_work, area_calc_chosen],
                }
                plot_vars_1d(oo_current)

    # - calculation of freq/rate at one stage -
    ws, gs = [], []
    if not flag_two_stages and flag_stat:
        for i_sample in range(n_samples):
            ids_one_t_interval, t_one_interval, _ = \
                mix.get_ids(t_work, t_intervals[i_sample])

            dict_wg = find_wg(
                t_one_interval,
                data_work[ids_one_t_interval],
                sel_wg,
                flag_print=False
            )
            w_est, line_w_est = give_res(dict_wg, 'est', 'w', coef_norm_global_w)
            g_est, line_g_est = give_res(dict_wg, 'est', 'g', coef_norm_global_g)
            w_adv, line_w_adv = give_res(dict_wg, 'adv', 'w', coef_norm_global_w)
            g_adv, line_g_adv = give_res(dict_wg, 'adv', 'g', coef_norm_global_g)

            if 'est' not in sel_wg:
                if w_adv is not None:
                    if np.abs((w_adv - w_est)/w_est) <= threshold_w:
                        ws.append(w_adv)
                if g_adv is not None:
                    if np.abs((g_adv - g_est) / g_est) <= threshold_g:
                        gs.append(g_adv)
            else:
                if w_est is not None:
                    ws.append(w_est)
                if g_est is not None:
                    gs.append(g_est)
        ws = np.array(ws)
        gs = np.array(gs)

    # - calculation of freq/rate at two stages -
    elif flag_two_stages and flag_stat:
        for i_sample in range(n_samples):

            # - FIND GAMMA -
            ids_one_t_interval, t_one_interval, _ = \
                mix.get_ids(t_gamma, t_intervals[i_sample])

            dict_gamma = find_wg(
                t_one_interval,
                data_gamma[ids_one_t_interval],
                sel_wg,
                flag_print=False
            )

            g_est, line_g_est = give_res(dict_gamma, 'est', 'g', coef_norm_global_g)
            g_adv, line_g_adv = give_res(dict_gamma, 'adv', 'g', coef_norm_global_g)
            if 'est' not in sel_wg:
                if g_adv is not None:
                    if np.abs((g_adv - g_est) / g_est) <= threshold_g:
                        gs.append(g_adv)
            else:
                if g_est is not None:
                    gs.append(g_est)

            # - FIND FREQUENCY -
            g_mult = g_est
            if g_adv is not None:
                g_mult = g_adv
            g_mult /= coef_norm_global_g
            dict_freq_filt = one_stage_filtering(
                t_work,  # filtering is performed in the work domain
                data_work * np.exp(- g_mult * t_work),
                oo_filt_freq
            )
            data_freq   = dict_freq_filt['filt']
            t_freq      = dict_freq_filt['x']

            ids_one_t_interval, t_one_interval, _ = \
                mix.get_ids(t_freq, t_intervals[i_sample])

            dict_freq = find_wg(
                t_one_interval,  # calculation is performed in one of the time domains
                data_freq[ids_one_t_interval],
                sel_wg,
                flag_print=False
            )

            w_est, line_w_est = give_res(dict_freq, 'est', 'w', coef_norm_global_w)
            w_adv, line_w_adv = give_res(dict_freq, 'adv', 'w', coef_norm_global_w)
            if 'est' not in sel_wg:
                if w_adv is not None:
                    if np.abs((w_adv - w_est)/w_est) <= threshold_w:
                        ws.append(w_adv)
            else:
                if w_est is not None:
                    ws.append(w_est)
        ws = np.array(ws)
        gs = np.array(gs)

    # - statistical results -
    stat_res = {
        'w': None, 'err_w': None,
        'g': None, 'err_g': None,
    }
    if flag_stat:
        # frequency histogram
        hist_w = np.histogram(ws, bins=n_bins)
        fit_mean_w, fit_sigma_w = stat_norm.fit(ws)
        f_data_w = stat_norm.pdf(hist_w[1], fit_mean_w, fit_sigma_w)

        # Rate histogram
        hist_g = np.histogram(gs, bins=n_bins)
        fit_mean_g, fit_sigma_g = stat_norm.fit(gs)
        f_data_g = stat_norm.pdf(hist_g[1], fit_mean_g, fit_sigma_g)

        err_w = GLO.COEF_ERR * fit_sigma_w
        err_g = GLO.COEF_ERR * fit_sigma_g

        # - plotting and printing statistical results -
        if flag_plot_print:
            # --- Frequency: Histogram ---
            # signal
            nsignals = 2  # histogram, normal distribution
            ch_signals = GLO.create_signals_dds(
                GLO.def_arbitrary_1d, [dd] * nsignals,
                flag_arbitrary=True,
                xs=[n_bins, hist_w[1]],
                datas=[ws, f_data_w],
            )

            # styling:
            ff_current = dict(ff)
            ff_current.update({
                'xlabel': '\omega[\omega_{ci}]',
                'ylabel': 'a.u',
                'title': 'Histogram:\ Frequency',
                'flag_maxlocator': True,
                'legends': ['histogram', 'normal\ distribution'],
                'flags_hist': [True, False],
            })

            # plotting
            oo_current = {
                'signals': ch_signals,
                'ff': ff_current,
                'flag_plot': flag_plot_internal,
            }
            curves_whist = plot_vars_1d(oo_current)

            # --- GAMMA: Histogram ---
            # signal
            nsignals = 2  # histogram, normal distribution
            ch_signals = GLO.create_signals_dds(
                GLO.def_arbitrary_1d, [dd] * nsignals,
                flag_arbitrary=True,
                xs=[n_bins, hist_g[1]],
                datas=[gs, f_data_g],
            )

            # styling:
            ff_current = dict(ff)
            ff_current.update({
                'xlabel': '\gamma[\omega_{ci}]',
                'ylabel': 'a.u',
                'title': 'Histogram:\ Damping/Growth Rate',
                'flag_maxlocator': True,
                'legends': ['histogram', 'normal\ distribution'],
                'flags_hist': [True, False],
            })

            # plotting
            oo_current = {
                'signals': ch_signals,
                'ff': ff_current,
                'flag_plot': flag_plot_internal,
            }
            curves_ghist = plot_vars_1d(oo_current)

            # plot histograms
            if flag_subplots:
                oo_sub = {
                    'ncols': 2,
                    'nrows': 1,
                    'list_curves': [curves_whist, curves_ghist],
                    'flag_subplots': flag_subplots,
                }
                plot_several_curves(oo_sub)

            # Print results
            line_stat = '--- STATISTICS ---\n'
            line_stat += 'number of frequency samples = ' + '{:d}'.format(len(ws)) + '\n'
            line_stat += 'number of rate samples = ' + '{:d}'.format(len(gs)) + '\n'
            line_stat += 'w' + line_norm_w + ' = ' + '{:0.3e}'.format(fit_mean_w) \
                        + '+-' + '{:0.3e}'.format(err_w) + '\n'
            line_stat += 'g' + line_norm_g + ' = ' + '{:0.3e}'.format(fit_mean_g) \
                        + '+-' + '{:0.3e}'.format(err_g)
            print(line_stat)

        # - save results of naive calculation -
        stat_res = {
            'w': fit_mean_w, 'err_w': err_w,
            'g': fit_mean_g, 'err_g': err_g,
        }
    out_res.update({'stat': stat_res})

    return out_res

