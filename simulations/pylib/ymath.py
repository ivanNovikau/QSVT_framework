import pylib.mix as mix
import pylib.Global_variables as GLO
import numpy as np
import scipy.signal
from scipy import stats
from scipy.optimize import curve_fit
from scipy import constants
from scipy import interpolate
import scipy.special
import sys


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(GLO)


def estimate_w(x, y):
    # find peaks of the signal
    ids_peaks, _ = scipy.signal.find_peaks(abs(y), height=0)
    if np.size(ids_peaks) < GLO.MIN_N_PEAKS:
        return None

    x_peaks = x[ids_peaks]
    y_peaks = abs(y[ids_peaks])

    # frequency
    w = np.pi / np.mean(np.diff(x_peaks))

    x_fit = x[ids_peaks[0]:ids_peaks[-1] + 1] - x[ids_peaks[0]]
    y_fit = y_peaks[0] * np.cos(w * x_fit)
    x_fit += x[ids_peaks[0]]

    # results
    out = {'ids_peaks': ids_peaks,
           'w': w,        'g': None,
           'w_err': None, 'g_err': None,
           'y_fit': y_fit, 'x_fit': x_fit}
    return out


def estimate_g(x, y):
    # initialize output:
    out = {}

    # find peaks of the signal
    ids_peaks, _ = scipy.signal.find_peaks(abs(y), height=0)
    if len(ids_peaks) >= GLO.MIN_N_PEAKS:
        x_peaks = x[ids_peaks]
        y_peaks = abs(y[ids_peaks])

        res = stats.theilslopes(
            np.log(y_peaks),
            x_peaks - x[ids_peaks[0]],
            GLO.CONFIDENCE_PERC
        )
        g = res[0]
        b0 = res[1]
        g_err_l = res[2]
        g_err_u = res[3]
        g_err = np.max([g_err_l, g_err_u])

        out['ids_peaks'] = ids_peaks

        x_fit = x[ids_peaks[0]:ids_peaks[-1] + 1] - x[ids_peaks[0]]
        y_fit = np.exp(b0 + g * x_fit)
        x_fit += x[ids_peaks[0]]
    else:
        out['ids_peaks'] = None

        res = stats.theilslopes(
            np.log(y),
            x - x[0],
            GLO.CONFIDENCE_PERC
        )
        g = res[0]
        b0 = res[1]
        g_err_l = res[2]
        g_err_u = res[3]
        g_err = np.max([g_err_l, g_err_u])

        x_fit = x - x[0]
        y_fit = np.exp(b0 + g * x_fit)
        x_fit += x[0]

    # results
    out['g'] = g
    out['w'] = None
    out['g_err'] = np.abs(g - g_err)
    out['w_err'] = None
    out['y_fit'] = y_fit
    out['x_fit'] = x_fit

    return out


def estimate_wg(x, y):
    # find peaks of the signal
    ids_peaks, _ = scipy.signal.find_peaks(abs(y), height=0)

    if len(ids_peaks) >= GLO.MIN_N_PEAKS:
        x_peaks = x[ids_peaks]
        y_peaks = abs(y[ids_peaks])

        # frequency
        w = np.pi / np.mean(np.diff(x_peaks))

        res = stats.theilslopes(
            np.log(y_peaks),
            x_peaks - x[ids_peaks[0]],
            GLO.CONFIDENCE_PERC
        )
        g = res[0]
        b0 = res[1]
        g_err_l = res[2]
        g_err_u = res[3]
        g_err = np.max([g_err_l, g_err_u])

        # fitting signals:
        x_fit = x[ids_peaks[0]:ids_peaks[-1] + 1] - x[ids_peaks[0]]
        y_fit = np.cos(w * x_fit) * np.exp(b0 + g * x_fit)
        x_fit += x[ids_peaks[0]]
        out = {'ids_peaks': ids_peaks,
               'w': w, 'g': g,
               'w_err': None, 'g_err': np.abs(g - g_err),
               'y_fit': y_fit, 'x_fit': x_fit}
    else:
        out = estimate_g(x, y)

    return out


def advanced_fitting(x, y, F, p0, id_w, id_g):
    popt, pcov = curve_fit(F, x, y, p0=p0)
    y_fit = F(x, *popt)
    wg_errs = np.sqrt(np.diag(pcov))

    if id_w is None:
        w, w_err = None, None
    else:
        id_w = np.int(id_w)
        w, w_err = popt[id_w], GLO.COEF_ERR*wg_errs[id_w]

    if id_g is None:
        g, g_err = None, None
    else:
        id_g = np.int(id_g)
        g, g_err = popt[id_g], GLO.COEF_ERR*wg_errs[id_g]

    res_adv = {
        'g': g,         'w': w,
        'g_err': g_err, 'w_err': w_err,
        'y_fit': y_fit, 'x_fit': x
    }
    return res_adv


def advanced_wg(x, y, flag_print=True):
    # estimation of w,g
    wg_est = estimate_wg(x, y)

    # initial result dictionary
    out = {
        'adv': None,
        'est': wg_est
    }

    try:
        if len(wg_est['ids_peaks']) < GLO.MIN_N_PEAKS:
            # fitting with only gamma
            F = lambda x_loc, A0, g: A0 * np.exp(g * x_loc)
            p0 = [y[0], wg_est['g']]
            out_adv = advanced_fitting(x, y, F, p0, None, 1)
        else:
            # fitting with gamma and frequency
            F = lambda x_loc, A0, g, w, Phase0: \
                A0 * np.cos(w * x_loc + Phase0) * np.exp(g * x_loc)
            p0 = [y[0], wg_est['g'], wg_est['w'], 0]
            out_adv = advanced_fitting(x, y, F, p0, 2, 1)
        out.update({'adv': out_adv})
    except:
        if flag_print:
            print('Advanced fitting with dynamic rate failed.')
        try:
            if len(wg_est['ids_peaks']) >= GLO.MIN_N_PEAKS:
                # fitting with frequency
                F = lambda x_loc, A0, w, Phase0: A0 * np.cos(w * x_loc + Phase0)
                p0 = [y[0], wg_est['w'], 0]
                out_adv = advanced_fitting(x, y, F, p0, 1, None)
                out.update({'adv': out_adv})
        except:
            if flag_print:
                print('Advanced fitting with frequency failed.')
    return out


def advanced_w(x, y, flag_print=True):
    # estimation of w,g
    wg_est = estimate_w(x, y)

    # initial result dictionary
    out = {
        'adv': None,
        'est': wg_est
    }

    # fitting with only frequency
    if len(wg_est['ids_peaks']) >= GLO.MIN_N_PEAKS:
        F = lambda x_loc, A0, w, Phase0: A0 * np.cos(w * x_loc + Phase0)
        p0 = [y[0], wg_est['w'], 0]
        try:
            out_adv = advanced_fitting(x, y, F, p0, 1, None)
            out.update({'adv': out_adv})
        except:
            if flag_print:
                print('Advanced fitting with frequency failed.')
    return out


def advanced_g(x, y, flag_print=True):
    # estimation of w,g
    wg_est = estimate_wg(x, y)

    # initial result dictionary
    out = {
        'adv': None,
        'est': wg_est
    }

    # fitting with only gamma
    F = lambda x_loc, A0, g: A0 * np.exp(g * x_loc)
    p0 = [y[0], wg_est['g'], None]

    try:
        out_adv = advanced_fitting(x, y, F, p0)
        out.update({'adv': out_adv})
    except:
        if flag_print:
            print('Advanced fitting with dynamic rate failed.')
    return out


def estimate_w_max_fft(w, f, oo={}):
    n_max = oo.get('n_max', 1)

    ids_peaks, _ = scipy.signal.find_peaks(f, height=0)
    f_peaks = f[ids_peaks]
    w_peaks = w[ids_peaks]
    ids_sort = np.argsort(f_peaks)

    w_max, f_max = np.zeros(n_max), np.zeros(n_max)
    nf = np.size(f_peaks)
    for count_max in range(n_max):
        id_sort = ids_sort[nf - 1 - count_max]
        w_max[count_max] = w_peaks[id_sort]
        f_max[count_max] = f_peaks[id_sort]

    out = {
        'w_max': w_max,
        'f_max': f_max
    }
    return out


def fft_y(x, y=None, oo=None):
    # for a 1d or 2d signal y

    if oo is None:
        oo = {}

    # additional parameters:
    flag_f2_arranged = oo.get('flag_f2_arranged', False)
    axis    = oo.get('axis', 0)  # 0, 1

    flag_yw = True
    if y is None:
        flag_yw = False

    # calculate frequency:
    nx2_ceil = np.ceil( np.log2(np.size(x)) )
    nx = int(2 ** nx2_ceil)
    nx_half = np.int(nx / 2)

    x_new = np.linspace(x[0], x[-1], nx)
    dx = np.min(np.diff(x_new))
    freq_max = 2*np.pi / dx
    dfreq = freq_max / nx

    w = np.array([dfreq * i for i in range(nx_half + 1)])

    left_a  = - np.flipud(w)
    right_a = w[1:np.size(w)-1]
    w2 = np.concatenate((left_a, right_a))

    # w2_ref = np.fft.fftfreq(nx, dx/(2*np.pi))  # to check in debuging the structure of the f2_raw

    # first part of results
    res = {
        'w': w, 'w2': w2, 'x_new': x_new
    }

    # calculate FFT of the signal y
    if not flag_yw:
        return res

    if nx != np.shape(y)[axis]:
        f_interp = interpolate.interp1d(x, y, axis=axis)
        y_new = f_interp(x_new)
    else:
        x_new, y_new = x, y
    f2_raw = np.fft.fft(y_new, n=nx, axis=axis)  # two-sided FFT
    f2 = np.abs(f2_raw / nx)

    f = None
    if np.size(np.shape(y)) == 1:
        f = f2[range(nx_half + 1)]
        f[1:np.size(f) - 1] = 2 * f[1:np.size(f) - 1]
    elif axis is 0:
        f = f2[0:nx_half + 1, :]
        f[1:np.size(f) - 1, :] = 2 * f[1:np.size(f) - 1, :]
    elif axis is 1:
        f = f2[:, 0:nx_half + 1]
        f[:, 1:np.size(f) - 1] = 2 * f[:, 1:np.size(f) - 1]

    if flag_f2_arranged:
        left_a = np.flipud(f2[0:nx_half+1])
        right_a = np.flipud(f2[nx_half+1:np.size(f2)])
        f2_arranged = np.concatenate((left_a, right_a))

        # update results
        res.update({'f2_arranged': f2_arranged})

    # update results:
    res.update({'f': f, 'f2_raw': f2_raw})

    return res

def filtering(x, y, oo):
    def res_None(x_loc, y_loc, norm_w_loc):
        ffres = fft_y(x_loc, y_loc, oo={'flag_f2_arranged': True})
        res = {'x': x,
               'filt': y_loc,
               'fft_init': ffres['f'],
               'fft_filt': None,
               'fft_init_2': ffres['f2_arranged'],
               'fft_filt_2': None,
               'w':  ffres['w']  * norm_w_loc,
               'w2': ffres['w2'] * norm_w_loc}
        return res

    # determine a filter type
    sel_filt = oo.get('sel_filt', None)  # None, 'rough', 'smooth', 'fft_smooth'
    norm_w   = oo.get('norm_w', 1)
    x_filt = np.array(x)
    y_filt = np.array(y)

    # no any filtering
    out = {}
    if sel_filt is None:
        out = res_None(x_filt, y_filt, norm_w)
    if sel_filt is 'rough':
        out = rough_filter(x_filt, y_filt, oo)
    if sel_filt is 'smooth':
        out = smooth_filter(x_filt, y_filt, oo)
    if sel_filt is 'fft_smooth':
        out = fft_smooth_filter(x_filt, y_filt, oo)

    out['x'] = x_filt

    return out


def rough_filter(x, y, oo):
    w_interval = oo.get('w_interval', None)
    norm_w = oo.get('norm_w', 1)

    ffres = fft_y(x, y, oo={'flag_f2_arranged': True})
    w,  w2  = ffres['w'], ffres['w2']

    w = w * norm_w
    w2 = w2 * norm_w

    nw, nw2 = np.size(w), np.size(w2)
    yw2 = ffres['f2_raw']

    # filtering
    if w_interval is not None:
        if w_interval[-1] is 0:
            yw2[0] = 0.0
        else:
            if w_interval[-1] is -1 or w_interval[-1] is None:
                w_interval[-1] = w[-1]

            ids_w, _, _ = mix.get_ids(w, w_interval)
            id_w1, id_w2 = ids_w[0], ids_w[-1]
            yw2[id_w1:id_w2 + 1] = 0.0
            yw2[nw2 - id_w2:nw2 - id_w1] = 0.0

    # new signal and its FFT
    y_new     = np.fft.irfft(yw2, np.size(ffres['x_new']))
    ffres_new = fft_y(ffres['x_new'], y_new, oo={'flag_f2_arranged': True})

    # build the new signal along the initial x-grid
    y_new = np.interp(x, ffres['x_new'], y_new)

    # results
    out = {'filt': y_new,
           'fft_init': ffres['f'],
           'fft_filt': ffres_new['f'],
           'fft_init_2': ffres['f2_arranged'],
           'fft_filt_2': ffres_new['f2_arranged'],
           'w': w, 'w2': w2}

    return out


def smooth_filter(x, y, oo):
    norm_w = oo.get('norm_w', 1)
    wind = int(oo.get('wind', 3))  # windows, which has to be an ODD number

    # smoothing
    y_filt = smooth(y, wind)

    # find FFT
    ffres      = fft_y(x, y,      oo={'flag_f2_arranged': True})
    ffres_filt = fft_y(x, y_filt, oo={'flag_f2_arranged': True})

    # results
    out = {'filt': y_filt,
           'fft_init': ffres['f'],
           'fft_filt': ffres_filt['f'],
           'fft_init_2': ffres['f2_arranged'],
           'fft_filt_2': ffres_filt['f2_arranged'],
           'w': ffres['w'] * norm_w,
           'w2': ffres['w2'] * norm_w
           }
    return out


def fft_smooth_filter(x, y, oo):
    norm_w = oo.get('norm_w', 1)
    w_interval = oo.get('w_interval', [0, 0])
    wind = oo.get('wind', 0)

    ffres = fft_y(x, y, oo={'flag_f2_arranged': True})
    w, w2   = ffres['w'], ffres['w2']
    nw, nw2 = np.size(w), np.size(w2)
    yw2     = ffres['f2_raw']

    w = w * norm_w
    w2 = w2 * norm_w

    # smoothing of the fft
    id_w1, _ = mix.find(w, w_interval[0])
    id_w2, _ = mix.find(w, w_interval[-1])
    yw2[id_w1:id_w2 + 1] = smooth(yw2[id_w1:id_w2 + 1], wind)
    yw2[nw2 - id_w2:nw2 - id_w1] = smooth(yw2[nw2 - id_w2:nw2 - id_w1], wind)

    # new signal and its FFT
    y_new     = np.fft.irfft(yw2, np.size(ffres['x_new']))
    ffres_new = fft_y(ffres['x_new'], y_new, oo={'flag_f2_arranged': True})

    # build the new signal along the initial x-grid
    if nw != np.size(y):
        y_new = np.interp(x, ffres['x_new'], y_new)

    # results
    out = {'filt': y_new,
           'fft_init': ffres['f'],
           'fft_filt': ffres_new['f'],
           'fft_init_2': ffres['f2_arranged'],
           'fft_filt_2': ffres_new['f2_arranged'],
           'w': w,
           'w2': w2
           }
    return out


def smooth(a, WSZ):
    # taken from
    # https://stackoverflow.com/questions/40443020/matlabs-smooth-implementation-n-point-moving-average-in-numpy-python

    # a: NumPy 1-D array containing the data to be smoothed
    # WSZ: smoothing window size needs, which must be ODD number,
    # as in the original MATLAB implementation
    out0 = np.convolve(a, np.ones(WSZ, dtype=int), 'valid') / WSZ
    r = np.arange(1, WSZ-1, 2)
    start = np.cumsum(a[:WSZ-1])[::2]/r
    stop = (np.cumsum(a[:-WSZ:-1])[::2]/r)[::-1]
    return np.concatenate((start, out0, stop))


def find_wc(B0, m, Z):
    # B0 in Tesla, m in kg, T in J
    wc = Z * constants.elementary_charge * B0 / m  # in rad/s
    return wc


def find_Te(Lx, a, B0, mf, Zf):
    # find temperature Ti of a species i:
    # tau_i = Ti / Te
    # mf - mass of the first species
    # Zf - charge (without e) of the first species
    # Remark: in ORB5 tau_e = 1 always
    wc = find_wc(B0, mf, Zf)
    Te = (2*a*wc/Lx)**2 * mf  # in J
    return Te


def find_Ti(tau_i, Lx, a, B0, mf, Zf):
    # find electron temperature Te:
    # mf - mass of the first species
    # Zf - charge (without e) of the first species
    # Remark: in ORB5 tau_e = 1 always
    Te = find_Te(Lx, a, B0, mf, Zf)
    Ti = Te * tau_i
    return Ti


def find_Lx(a, Te, B0, mf, Zf):
    # mf - mass of the first species (kg)
    # Zf - charge (without e) of the first species
    # Te - electron temperature (J)
    # a - minor radius (m)
    # B0 - magnetic field (T)
    wc = find_wc(B0, mf, Zf)
    cs = find_cs(Te, mf)
    Lx = 2*a*wc / cs
    return Lx


def find_vt(T, m):
    # m in kg, T in J
    vt = np.sqrt(T/m)
    return vt


def find_cs(Te, m):
    # m in kg, Te in J
    cs = np.sqrt(Te/m)
    return cs


def find_k(krpert, Tf, B0, Lwork, mf, Zf):
    # Find normalised radial wavenumber
    # mf - mass of the first species (kg)
    # Zf - charge (without e) of the first species
    # Tf - temperature of the first species (J)
    # B0 - magnetic field (Tesla)
    # krpert - radial vector from an input file
    rhoL = find_rhoL(Tf, B0, mf, Zf)
    k = krpert * (np.pi / Lwork) * rhoL
    return k


def find_rhoL(T, B0, m, Z):
    # Find a value of the Larmor radius
    # m - mass of a species (kg)
    # Z - charge (without e) of a species
    # T - temperature of a species (J)
    # B0 - magnetic field (Tesla)
    wc = find_wc(B0, m, Z)
    vt = find_vt(T, m)
    rhoLarmor = np.sqrt(2) * vt / wc
    return rhoLarmor


def find_rho_star(Lx):
    return 2. / Lx


def find_norm(y, y_norm_to=None):
    if y_norm_to is None:
        abs_y = np.abs(y)
    else:
        abs_y = np.abs(y_norm_to)
    y_norm = y / np.nanmax(abs_y[abs_y != np.inf])
    return y_norm


def get_v_res(dd, s1, w_wci, sel_species):
    cs_speak = np.sqrt(dd['electrons'].T_speak(dd) / dd['pf'].mass)
    norm_v = np.max(dd[sel_species].nT_equil['T']) / cs_speak

    # Te_max = np.max(dd['electrons'].nT_equil['T_J'])
    # cs_max = np.sqrt(Te_max / dd['pf'].mass)
    # norm_v = 1 / cs_max

    w0 = w_wci * dd['wc']

    q_s1 = mix.get_data_at_x1(s1, dd['q']['s'], dd['q']['data'])
    vres = (q_s1 * dd['R0'] * w0) * norm_v

    return vres


def get_w_from_v_res(dd, s1, vres, sel_species):
    cs_speak = np.sqrt(dd['electrons'].T_speak(dd) / dd['pf'].mass)
    norm_v = np.max(dd[sel_species].nT_equil['T']) / cs_speak

    q_s1 = mix.get_data_at_x1(s1, dd['q']['s'], dd['q']['data'])
    w0 = vres / (q_s1 * dd['R0'] * norm_v)
    w_wci = w0 / dd['wc']

    return w_wci


def avr_x1x2(vvar, one_signal, oo ):
    # one_signal = {'avr_operation': sel_av, 'avr_domain': av_domain}
    # sel_av = 'type_av-coord_av'
    # type_av = 'mean', 'rms', 'point'
    # coord_av = 't', 's', 'vpar', 'mu', 'chi', 'n', 'm'
    # av_domain is x_point or x_interval

    # type of averaging and the coordinate along which the averaging will be performed
    sel_av = one_signal['avr_operation']
    type_av, coord_av = sel_av.split('-')
    if type_av == 'none':
        vvar_avr = dict(vvar)
        vvar_avr['line_avr'] = ''
        return vvar_avr
    flag_aver_interp = oo.get('flag_aver_interp', None)

    # data and coordinates
    data = vvar['data']
    x1 = vvar[vvar['x1']]
    x2 = vvar[vvar['x2']]

    # define a coordinate axis to average along:
    if coord_av == vvar['x1']:
        dir_av = 0
        x_av   = x1
        format_x_av = vvar['fx1']
        x_work = x2
        format_x_work = vvar['fx2']
        name_x_work = vvar['x2']
    elif coord_av == vvar['x2']:
        dir_av = 1
        x_av   = x2
        format_x_av = vvar['fx2']
        x_work = x1
        format_x_work = vvar['fx1']
        name_x_work = vvar['x1']
    else:
        mix.error_mes('Wrong name of a coordinate for averaging.')

    # use reference axis for averaging:
    if flag_aver_interp:
        x_av_ref = oo.get(coord_av + '_ref', None)
        f_interp = interpolate.interp1d(x_av, data, axis=dir_av, fill_value="extrapolate")
        data = f_interp(x_av_ref)
        x_av = np.array(x_av_ref)

    # averaging domain
    av_domain = one_signal.get('avr_domain', [x_av[0], x_av[-1]])

    # --- averaging ---
    data_av, line_av, opt_av = None, '', None
    ids_av, _, line_x = mix.get_ids(x_av, av_domain, format_x=format_x_av)
    if coord_av == vvar['x1']:
        temp = data[ids_av, :]
    else:
        temp = data[:, ids_av]

    if type_av == 'mean':
        data_av = np.mean(temp, axis=dir_av)
        line_av = 'mean_' + coord_av + ':\ '
        line_av += coord_av + ' = ' + line_x
    elif type_av == 'rms':
        data_av = np.sqrt(np.mean(temp ** 2, axis=dir_av))
        line_av = 'rms_' + coord_av + ':\ '
        line_av += coord_av + ' = ' + line_x
    elif type_av == 'point':
        data_av = np.squeeze(temp)
        line_av += coord_av + ' = ' + line_x
    elif type_av == 'sum':
        data_av = np.sum(temp, axis=dir_av)
        line_av = 'sum_' + coord_av + ':\ '
        line_av += coord_av + ' = ' + line_x
    elif type_av == 'max':
        data_av = np.max(temp, axis=dir_av)
        opt_av = x_av[np.argmax(temp, axis=dir_av)]
        line_av = 'max_{' + coord_av + '}'
    elif type_av == 'absmax':
        abs_temp  = np.abs(temp)
        data_av = np.max(abs_temp, axis=dir_av)
        opt_av  = x_av[np.argmax(abs_temp, axis=dir_av)]
        line_av = 'max_{' + coord_av + '}'
    else:
        mix.error_mes('Wrong name of an averaging operation.')

    # form result dictionary
    vvar_avr = {
        'data':        data_av,
        'x':           x_work,
        'name_x':      name_x_work,
        'format_x':    format_x_work,
        'line_avr':    line_av,
        'opt_av':      opt_av
    }

    # additional fields:
    name_fields = ['tit', 'labx', 'laby']
    for name_field in name_fields:
        if name_field in vvar:
            vvar_avr[name_field] = vvar[name_field]

    return vvar_avr


def post_processing(data, x, oo_operations):
    # - check operations -
    if oo_operations is None:
        return data, x

    # -- PERFORM CONSEQUENTLY ALL OPERATIONS --
    data_work, x_work = np.array(data), np.array(x)
    for oo_operation in oo_operations:
        # - operation type and working domain -
        sel_operation = oo_operation.get('operation', None)
        x_domain      = oo_operation.get('domain', None)
        if x_domain is None:
            x_domain = [x_work[0], x_work[-1]]
        data_work, x_work = mix.get_x_data_interval(x_domain, x_work, data_work)
        data_temp, x_temp = [], []

        # * No operations *
        if sel_operation is None:
            data_temp, x_temp = data_work, x_work

        # * integration from beginning till the point for every selected x-steps *
        elif sel_operation == 'integration-accumulation':
            width_domain = oo_operation.get('width', None)
            if width_domain is None:
                mix.error_mes('Postprocessing: You need not-None \'width\' field '
                              'for \'integration-accumulation\' postprocessing')

            x_start = x_work[0]
            x_end = x_start + width_domain
            while x_end <= x_work[-1]:
                data_curr, x_curr = mix.get_x_data_interval([x_start, x_end], x_work, data_work)
                value_curr = np.trapz(data_curr, x=x_curr)

                data_temp.append(value_curr)
                x_temp.append(x_end)

                x_end += width_domain

        # * integration within every x-step *
        elif sel_operation == 'integration-step':
            width_domain = oo_operation.get('width', None)
            step_domain  = oo_operation.get('step', None)
            if width_domain is None:
                mix.error_mes('Postprocessing: You need not-None \'width\' field '
                              'for \'integration-step\' postprocessing')
            if step_domain is None:
                mix.error_mes('Postprocessing: You need not-None \'step\' field '
                              'for \'integration-step\' postprocessing')

            x_start = x_work[0]
            x_end   = x_start + width_domain
            while x_end <= x_work[-1]:
                data_curr, x_curr = mix.get_x_data_interval([x_start, x_end], x_work, data_work)
                value_curr = np.trapz(data_curr, x=x_curr)

                x_point = (x_end + x_start) / 2.
                data_temp.append(value_curr)
                x_temp.append(x_point)

                x_start += step_domain
                x_end    = x_start + width_domain

        # * get peaks of absolute signal *
        elif sel_operation == 'abs_peaks':
            ids_abs_peaks, _ = scipy.signal.find_peaks(np.abs(data_work))
            x_temp    = x_work[ids_abs_peaks]
            data_temp = data_work[ids_abs_peaks]

        # * take absolute signal *
        elif sel_operation == 'abs':
            data_temp, x_temp = np.abs(data_work), x_work

        # * multiply by a value *
        elif sel_operation == 'mult':
            coef = oo_operation.get('coef', None)
            if coef is None:
                mix.error_mes('Postprocessing: You need not-None \'coef\' field '
                              'for \'mult\' postprocessing')

            data_temp, x_temp = coef * data_work, x_work

        # * multiply x by a value *
        elif sel_operation == 'mult-x':
            coef = oo_operation.get('coef', None)
            if coef is None:
                mix.error_mes('Postprocessing: You need not-None \'coef\' field '
                              'for \'mult\' postprocessing')

            data_temp, x_temp = data_work, coef * x_work

        # * filtering *
        elif sel_operation == 'filtering':
            oo_filters = oo_operation.get('oo_filters', [GLO.NONE_FILTER])
            data_filt, x_filt = global_filter_signal(x_work, data_work, oo_filters)
            x_temp = x_work
            data_temp = np.interp(x_temp, x_filt, data_filt)

        # * 1d Fourier transformation *
        elif sel_operation == 'fft-1d':
            oo_fft = oo_operation.get('oo_fft', [GLO.DEF_FFT])
            data_fft, x_fft = get_fft_1d(x_work, data_work, oo_fft)
            data_temp, x_temp = data_fft, x_fft

        # * Shift along the x *
        elif sel_operation == 'shift-x':
            x_shift_start = oo_operation.get('x-shift-start', None)
            id_x_shift_start, x_shift_start, _ = \
                mix.get_ids(x_work, x_shift_start)
            x_temp = x_work[id_x_shift_start:-1] - x_shift_start
            data_temp = data_work[id_x_shift_start:-1]

        # * wrong operation name *
        else:
            mix.error_mes('Wrong operation name')

        # - save modification -
        data_work, x_work = np.array(data_temp), np.array(x_temp)

    return data_work, x_work


def post_processing_2d(data, x, y, name_x, name_y, oo_operations):
    # - check operations -
    if oo_operations is None:
        return data, x, name_x, y, name_y

    # names of variables:
    name_x_res, name_y_res = name_x, name_y

    # -- PERFORM CONSEQUENTLY ALL OPERATIONS --
    data_work, x_work, y_work = np.array(data), np.array(x), np.array(y)
    for oo_operation in oo_operations:
        # - operation type and working domain -
        sel_operation = oo_operation.get('operation', None)
        x_domain = oo_operation.get('domain-' + name_x_res, None)
        y_domain = oo_operation.get('domain-' + name_y_res, None)
        if x_domain is None:
            x_domain = [x_work[0], x_work[-1]]
        if y_domain is None:
            y_domain = [y_work[0], y_work[-1]]
        data_work, x_work, y_work = mix.get_data_interval_2d(
            data, x_work, x_domain, y_work, y_domain
         )
        data_temp, x_temp, y_temp = [], [], []

        # * No operations *
        if sel_operation is None:
            data_temp, x_temp, y_temp = data_work, x_work, y_work

        # * multiply the function by a value *
        elif sel_operation == 'mult':
            coef = oo_operation.get('coef', None)
            if coef is None:
                mix.error_mes('Postprocessing: You need not-None \'coef\' field '
                              'for \'mult\' postprocessing')
            data_temp, x_temp, y_temp = coef * data_work, x_work, y_work

        # * multiply the axis X by a value *
        elif sel_operation == 'mult-axis-x':
            coef = oo_operation.get('coef', None)
            if coef is None:
                mix.error_mes('Postprocessing: You need not-None \'coef\' field '
                              'for \'mult-axis-x\' postprocessing')
            data_temp, x_temp, y_temp = data_work, coef * x_work, y_work

        # * take an absolute signal *
        elif sel_operation == 'abs':
            data_temp, x_temp, y_temp = np.abs(data_work), x_work, y_work

        # * 2d Fourier transformation *
        elif sel_operation == 'fft-2d':
            oo_fft = oo_operation.get('oo_fft', [GLO.DEF_FFT_2D])
            data_work, x_work, name_x_res, y_work, name_y_res = get_fft_2d(
                x_work, name_x_res, y_work, name_y_res, data_work, oo_fft
            )
            data_temp, x_temp, y_temp = data_work, x_work, y_work

        # * wrong operation name *
        else:
            mix.error_mes('Wrong operation name')

        # - save modification -
        data_work, x_work, y_work = \
            np.array(data_temp), np.array(x_temp), np.array(y_temp)

    return data_work, x_work, name_x_res, y_work, name_y_res


def global_several_filters(x, y, oo_filt_loc):
    oo_filts_res = []
    if oo_filt_loc is None or len(oo_filt_loc) is 0:
        oo_filts_res.append(GLO.NONE_FILTER)
    elif isinstance(oo_filt_loc, dict):
        oo_filts_res.append(oo_filt_loc)
    else:
        oo_filts_res = oo_filt_loc  # array of filters

    # apply one by one the filters in the array :
    dict_loc = {'x': np.array(x), 'filt': np.array(y)}
    for one_filt in oo_filts_res:
        dict_loc = filtering(dict_loc['x'], dict_loc['filt'], one_filt)
    return dict_loc


def global_filter_signal(t, y, oo_filt):
    # initial fft
    init_dict = {'x': t, 'data': y}

    # filtering
    res_dict = dict(init_dict)
    filt = global_several_filters(init_dict['x'], init_dict['data'], oo_filt)
    res_dict['data'] = np.array(filt['filt'])
    res_dict['x']    = np.array(filt['x'])

    return res_dict['data'], res_dict['x']


def get_fft_1d(x, data, oo):
    # one or two-sided FFT
    flag_f2 = oo['flag_f2']

    # - frequency grid -
    ffres = fft_y(x)
    w = ffres['w2'] if flag_f2 else ffres['w']

    # - fft -
    var_fft = fft_y(x, data, {'flag_f2_arranged': flag_f2})
    data_fft = var_fft['f2_arranged'] if flag_f2 \
        else var_fft['f']

    return data_fft, w


def get_fft_2d(x, name_x, y, name_y, data, oo):
    # compare name_x and name_y with oo.name_coord_fft
    # if name_x = name_coord_fft -> transpose the data to gurantee
    #   that the FFT will be perform along the second axis

    # name of the coordinate axis, along which FFT will be taken
    name_coord_fft = oo.get('name_coord_fft', None)

    # one or two-sided FFT
    flag_f2 = oo.get('flag_f2', False)

    # axis along which will be performed the FFT
    axis_fft = 1

    # name of a new axis where FFT has been performed
    name_w = 'w'

    # choose coordinate to perform Fourier transformation
    if name_coord_fft == name_x:
        coord_fft = x
        coord_along = y
        data = data.T
        name_along = name_y
    else:
        coord_fft = y
        coord_along = x
        name_along = name_x

    # - frequency grid -
    ffres = fft_y(coord_fft)
    w = ffres['w2'] if flag_f2 else ffres['w']

    # - FFT -
    data_fft = fft_y(coord_fft, data,
                    {'flag_f2_arranged': flag_f2, 'axis': axis_fft})

    # - result FFT -
    data_fft = data_fft['f2_arranged'] if flag_f2 else data_fft['f']

    return data_fft, coord_along, name_along, w, name_w


def decrease_mesh_size(X, Y, Z, step_x, step_y):
    print('Initial number of rows in X,Y,Z is: {:d}, {:d}, {:d}'.format(
        len(X), len(Y), len(Z))
    )
    print('Initial number of columns in X,Y,Z is: {:d}, {:d}, {:d}'.format(
        len(X[0]), len(Y[0]), len(Z[0]))
    )
    print('Initial size is {:d}'.format(np.size(Z)))

    X_red = np.array(X[::step_x, ::step_y])
    Y_red = np.array(Y[::step_x, ::step_y])
    Z_red = np.array(Z[::step_x, ::step_y])

    print('Result size is {:d}'.format(np.size(Z_red)))
    print('Result number of rows in Z is {:d}'.format( len(Z_red)) )
    print('Result number of columns in Z is {:d}'.format( len(Z_red[0])) )

    return X_red, Y_red, Z_red


def build_Wigner(dd):
    print("Computing... ", end='\r')

    data = dd['data']
    x_grid_orig = dd['x-grid-orig']
    k_grid = dd['k-grid']
    ids_work = dd['ids-x-work']
    x_ref = dd['x-ref']
    k_ref = dd['k-ref']

    Nx = len(ids_work)
    Nk = len(k_grid)

    # use this minimum as a maximum window for the Wigner function
    x_step   = x_grid_orig[1] - x_grid_orig[0]
    N_ids_dist = np.min([ids_work[0], len(x_grid_orig) - ids_work[-1]])

    ww = np.zeros([Nx, Nk], dtype = np.complex)
    for ix in range(len(ids_work)):
        sys.stdout.flush()
        if ix%20 == 0:
            print("Computing... id-x = {:d}".format(ix), end='\r')
        id_x1 = ids_work[ix]
        for ik in range(len(k_grid)):
            k1 = k_grid[ik]
            ww_k1 = 0
            for id_dist in range(N_ids_dist):
                dist = 2 * x_step * id_dist
                # ww_k1 += np.exp(-1j*k1*dist) * data[id_x1 + id_dist] * np.conj(data[id_x1 - id_dist])
                ww_k1 += np.exp(1j*k1*dist) * data[id_x1 + id_dist] * np.conj(data[id_x1 - id_dist])
            ww_k1 = 1./(2*np.pi*N_ids_dist) * ww_k1
            ww[ix, ik] = ww_k1

    res = {
        'data': ww,
        'ids-x': ids_work,
        'k': k_grid
    }

    print()
    print("Done.")

    return res










