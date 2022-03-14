from warnings import resetwarnings
import numpy as np
import scipy.special as ssp
import scipy.constants as sc
import matplotlib.pyplot as plt
from ipywidgets import interact, widgets
import scipy.interpolate 
import pylib.mix as mix
plt.rcParams.update({
    "text.usetex": True,
    # "font.family": "sans-serif",
    # "font.sans-serif": ["Helvetica"],
    'text.latex.preamble': r"\usepackage{amsmath} \boldmath"
})

def reload():
    mix.reload_module(mix)
    
# constants in Gauss units
me = sc.m_e * 1e3 # in gramms
qe = sc.e * 3e9  # in statcoulomb
c_light = sc.c * 1e2 # cm/s
pi4 = 4*np.pi

def print_constants_Gauss():
    print("electron mass: {:0.3e} g".format(me))
    print("electron charge: {:0.3e} statcoulomb".format(qe))
    print("light speed: {:0.3e} cm/s".format(c_light))
    return

def wp_func(n, q=qe, m = me):
    return np.sqrt(4*np.pi * n * q**2 / m)

def wc_func(B, q=qe, m = me):
    return q*B/(m*c_light)

def wuh_func(n, B, q=qe, m = me):
    return np.sqrt(wp_func(n,q,m)**2 + wc_func(B,q,m)**2)

def wr_func(n, B, q=qe, m = me):
    wc1 = wc_func(B,q,m)
    wp1 = wp_func(n,q,m)
    return (- wc1 + np.sqrt(wc1**2 + 4*wp1**2))/2

def wl_func(n, B, q=qe, m = me):
    wc1 = wc_func(B,q,m)
    wp1 = wp_func(n,q,m)
    return (wc1 + np.sqrt(wc1**2 + 4*wp1**2))/2

def dispersion_equ_X(w, kx_input, n, B, q=qe, m = me):
    wc1, wp1 = wc_func(B,q,m), wp_func(n,q,m)
    temp1 = kx_input**2 * c_light**2 * (w**2 - wuh_func(n, B)**2)
    temp2 = w**2 + wc1 * w - wp1**2
    temp3 = w**2 - wc1 * w - wp1**2
    return temp1 - temp2*temp3

def profile_n(n_center, x, x_n0, x_nw, n_center_2 = np.nan, x_n0_2 = np.nan, x_nw_2 = np.nan):
    y = (x - x_n0) / x_nw

    # --- option 1 ---
    n = n_center * np.exp(- y**2 / 2)
    if not np.isnan(n_center_2) and not np.isnan(x_n0_2) and not np.isnan(x_nw_2):
        y2 = (x - x_n0_2) / x_nw_2
        n2 = n_center_2 * np.exp(- y2**2 / 2)
        n += n2

    # --- option 2 ---
    # n = n_center * np.exp(- y**2 / 2) * (1 + ssp.erf(-0.3*y))
    # n0 = np.max(n)
    # n = n/n0 * n_center

    # --- option 3 ---
    # n = np.array([n_center]*len(x))

    # --- option 4 ---
    # n = n_center*(1.02 + np.cos(np.pi*x_prof/r0))/2 # DIII
    
    return n

def profile_B(B0, R0, x, B0_2=np.nan, R0_2=np.nan, x_B_left=np.nan, n_gap=10):
    B = B0 * R0 / (R0 + x)

    if not np.isnan(B0_2) and not np.isnan(R0_2) and not np.isnan(x_B_left):
        id_x_left, x_B_left = mix.find(x, x_B_left)

        x_left =  x[0:(id_x_left - n_gap)] 
        B2_left = np.array([0] * len(x_left))

        x_right = x[id_x_left:len(x)]
        B2_right = B0_2 * R0_2 / (R0_2 + (x_right - x_B_left))

        x2 = np.concatenate((x_left, x_right))
        B2 = np.concatenate((B2_left, B2_right))

        f_interp = scipy.interpolate.interp1d(x2, B2, kind='cubic')

        B2 = np.array(f_interp(x))

        B += B2

    return B

def normalization(n0, B0, sel_norm):
    va = B0 / np.sqrt(pi4*n0*me)
    wp = wp_func(n0)

    if sel_norm == 'va':
        coef_v_norm = va
        coef_t_norm = 1/wp
        coef_n_norm = n0
    elif sel_norm == 'c':
        coef_v_norm = c_light
        coef_t_norm = 1/wp
        coef_n_norm = n0
    else:
        print("wrong normalization selector")
        exit(-1)
    
    coef_r_norm = coef_v_norm * coef_t_norm
    coef_e_norm = coef_r_norm * pi4 * n0 * qe 

    return coef_v_norm, coef_t_norm, coef_r_norm, coef_e_norm, coef_n_norm

# 1-D: density and background magnetic field profiles + normalization + boundary profiles:
def prepare_norm_profiles_outgoing(dd, flag_full_print=False):
    # input parameters:
    r0 = dd["r0"]
    t0 = dd["t0"]
    n0 = dd["n0"]
    n0_2 = dd["n0_2"]
    B0 = dd["B0"]
    R0 = dd["R0"]
    B0_2 = dd["B0_2"]
    R0_2 = dd["R0_2"]
    x_B_left = dd["x_B_left"]
    kx = dd["kx"]
    nx = dd["nx"]
    nt = dd["nt"]
    x_n0 = dd["x_n0"]
    x_nw = dd["x_nw"]
    x_n0_2 = dd["x_n0_2"]
    x_nw_2 = dd["x_nw_2"]
    sel_norm = dd["sel-norm"]
    coef_coupling = dd["coef-coupling"]
    np_gap_B = dd["np-gap-B"]

    # coefficient to take into account additional multipliers in an oracle due to superpositions:
    coef_superposition = dd["coef-superposition"] 

    # normalization factors:
    coef_v_norm, coef_t_norm, coef_r_norm, coef_e_norm, coef_n_norm = normalization(n0, B0, sel_norm)

    # space coordinate:
    x = np.linspace(-r0, r0, nx) 
    x_norm = x / coef_r_norm
    x_plot = x / r0
    h = x_norm[1] - x_norm[0]
    dx = h * coef_r_norm  # in cm

    wavelen = 2*np.pi/kx

    # time interval: - tau, 0, tau, 2*tau, ...
    tau = 1. * t0 / (nt - 2)   # dt_norm
    t_norm = [(ii - 1)*tau for ii in range(nt)]
    dt = tau * coef_t_norm  # in seconds

    # reference frequency:
    w_ref  = kx * c_light  # rad/s
    w_ref_norm = w_ref * coef_t_norm
    
    # profiles
    n = profile_n(n0, x, x_n0, x_nw, n0_2, x_n0_2, x_nw_2)
    B = profile_B(B0, R0, x, B0_2, R0_2, x_B_left, np_gap_B)

    # normalized variables:
    n_norm = n / coef_n_norm
    c_norm = c_light / coef_v_norm
    B_norm = B / coef_e_norm
    kx_norm = kx * coef_r_norm
    kx_plot = kx * r0

    B0_norm = B0/coef_e_norm

    # Courant number 
    C = c_norm*tau/h

    # Hamiltonian normalization
    Bmax_norm = np.max(B_norm)
    nmax_norm = np.max(n_norm)

    # QC normalization of a Hamiltonian
    # norm_of_h = np.sqrt(Bmax_norm**2 + nmax_norm + 1./(2*h**2) + 
    #     coef_coupling**2 + w_ref_norm**2) * coef_superposition**2

    norm_of_h = np.sqrt(Bmax_norm**2 + nmax_norm + 1./(2*h**2) + 
        coef_coupling**2 + w_ref_norm**2) * coef_superposition**2

    # print parameters:
    if flag_full_print:
        print("x-step: {:0.3e} mm".format(dx*10))
        print("wavelength: {:0.3e} mm".format(wavelen * 10))
        print("wavelength/r0: {:0.3e}".format(wavelen/r0))
        print("n of xsteps in one wavelength: {:d}".format(int(wavelen/dx)))
        print("")
        print("time step: {:0.3e} s".format(dt))
        print("T-ref = 2pi/(kc): {:0.3e} s".format(2*np.pi/w_ref))
        print("T-ref-norm: {:0.3e}".format(2*np.pi/(w_ref*coef_t_norm)))
        print("n of tsteps in one period: {:d}".format(int(2*np.pi/(w_ref*dt))) )
        print("norm w-ref: {:0.3e}".format(w_ref_norm))
        print("\n------------------------------")
        print("--- Other variable ---")
        print("c-light norm: {:0.3e}".format(c_norm))
        print("")
        print("coef-v-norm: {:0.3e} cm/s".format(coef_v_norm))
        print("coef-t-norm: {:0.3e} s".format(coef_t_norm))
        print("coef-r-norm: {:0.3e} cm".format(coef_r_norm))
        print("coef-e-norm: {:0.3e} statvolt/cm or gauss".format(coef_e_norm))
        print("coef-n-norm: {:0.3e} cm-3".format(coef_n_norm))
        if sel_norm == "c":
            print("")
            print("c * sqrt(4pi nm): {:0.3e}".format(c_light * np.sqrt(pi4*n0*me)))
    
    # results:
    res = {
        'x': x, 'x-norm': x_norm, 'x-plot': x_plot, 'dx': dx, 'dx-norm': h,
        't-norm': t_norm, 'dt-norm': tau,
        'n': n, 'B': B, 'n-norm': n_norm, 'B-norm': B_norm,
        'B0-norm': B0_norm,
        'c-norm': c_norm,
        'kx-norm': kx_norm, 'kx-plot': kx_plot,
        'w-ref':   w_ref,   'w-ref-norm': w_ref_norm,
        'coef-v-norm': coef_v_norm, 'coef-t-norm': coef_t_norm, 'coef-r-norm': coef_r_norm,  
        'coef-e-norm': coef_e_norm, 'coef-n-norm': coef_n_norm,
        'norm-of-h': norm_of_h, "1/2h": 1./(2*h),
        "Courant": C
    }

    res.update(dd)

    return res

# 1-D: density and background magnetic field profiles + normalization:
def prepare_system(dd, flag_full_print=False):
    # input parameters:
    r0 = dd["r0"]
    t0 = dd["t0"]
    n0 = dd["n0"]
    B0 = dd["B0"]
    R0 = dd["R0"]
    kx = dd["kx"]
    nx = dd["nx"]
    nt = dd["nt"]
    x_n0 = dd["x_n0"]
    x_nw = dd["x_nw"]
    sel_norm = dd["sel-norm"]
    coef_coupling = dd["coef-coupling"]

    # number of hamiltonians to get necessary superpositions:
    coef_superposition = dd["coef-superposition"] 

    # normalization factors:
    coef_v_norm, coef_t_norm, coef_r_norm, coef_e_norm, coef_n_norm = normalization(n0, B0, sel_norm)

    # space coordinate:
    x = np.linspace(-r0, r0, nx) 
    x_norm = x / coef_r_norm
    x_plot = x / r0
    h = x_norm[1] - x_norm[0]
    dx = h * coef_r_norm  # in cm

    wavelen = 2*np.pi/kx

    # time interval: - tau, 0, tau, 2*tau, ...
    tau = 1. * t0 / (nt - 2)   # dt_norm
    t_norm = [(ii - 1)*tau for ii in range(nt)]
    dt = tau * coef_t_norm  # in seconds

    # reference frequency:
    w_ref  = kx * c_light  # rad/s
    w_ref_norm = w_ref * coef_t_norm
    
    # profiles
    n = profile_n(n0, x, x_n0, x_nw)
    B = profile_B(B0, R0, x)

    # normalized variables:
    n_norm = n / coef_n_norm
    c_norm = c_light / coef_v_norm
    B_norm = B / coef_e_norm
    kx_norm = kx * coef_r_norm
    kx_plot = kx * r0

    B0_norm = B0/coef_e_norm

    # Courant number 
    C = c_norm*tau/h

    # Hamiltonian normalization
    Bmax_norm = np.max(B_norm)
    nmax_norm = np.max(n_norm)

    # QC normalization of a Hamiltonian
    # norm_of_h = np.sqrt(Bmax_norm**2 + nmax_norm + 1./(2*h**2) + 
    #     coef_coupling**2 + w_ref_norm**2) * coef_superposition**2

    norm_of_h = np.sqrt(Bmax_norm**2 + nmax_norm + 1./(2*h**2) + 
        coef_coupling**2 + w_ref_norm**2) * coef_superposition

    # print parameters:
    if flag_full_print:
        print("x-step: {:0.3e} mm".format(dx*10))
        print("wavelength: {:0.3e} mm".format(wavelen * 10))
        print("wavelength/r0: {:0.3e}".format(wavelen/r0))
        print("n of xsteps in one wavelength: {:d}".format(int(wavelen/dx)))
        print("")
        print("time step: {:0.3e} s".format(dt))
        print("T-ref = 2pi/(kc): {:0.3e} s".format(2*np.pi/w_ref))
        print("T-ref-norm: {:0.3e}".format(2*np.pi/(w_ref*coef_t_norm)))
        print("n of tsteps in one period: {:d}".format(int(2*np.pi/(w_ref*dt))) )
        print("norm w-ref: {:0.3e}".format(w_ref_norm))
        print("\n------------------------------")
        print("--- Other variable ---")
        print("c-light norm: {:0.3e}".format(c_norm))
        print("")
        print("coef-v-norm: {:0.3e} cm/s".format(coef_v_norm))
        print("coef-t-norm: {:0.3e} s".format(coef_t_norm))
        print("coef-r-norm: {:0.3e} cm".format(coef_r_norm))
        print("coef-e-norm: {:0.3e} statvolt/cm or gauss".format(coef_e_norm))
        print("coef-n-norm: {:0.3e} cm-3".format(coef_n_norm))
        if sel_norm == "c":
            print("")
            print("c * sqrt(4pi nm): {:0.3e}".format(c_light * np.sqrt(pi4*n0*me)))
    
    # results:
    res = {
        'x': x, 'x-norm': x_norm, 'x-plot': x_plot, 'dx': dx, 'dx-norm': h,
        't-norm': t_norm, 'dt-norm': tau,
        'n': n, 'B': B, 'n-norm': n_norm, 'B-norm': B_norm,
        'B0-norm': B0_norm,
        'c-norm': c_norm,
        'kx-norm': kx_norm, 'kx-plot': kx_plot,
        'w-ref':   w_ref,   'w-ref-norm': w_ref_norm,
        'coef-v-norm': coef_v_norm, 'coef-t-norm': coef_t_norm, 'coef-r-norm': coef_r_norm,  
        'coef-e-norm': coef_e_norm, 'coef-n-norm': coef_n_norm,
        'norm-of-h': norm_of_h, "1/2h": 1./(2*h),
        "Courant": C
    }

    res.update(dd)

    return res


def plot_nB(ax, x, n, B, flag_norm):
    label_x, label_n, label_b = 'x (cm)', '$n[cm^{-3}]$', '$B[Gauss]$'
    title = "Profiles in real units"
    name_line_n, name_line_B = 'n', 'B'
    if flag_norm:
        label_x, label_n, label_b = 'x/r0', '$n/n_0$', '$B/B_0$'
        title = "Normalized profiles"
        name_line_n, name_line_B = 'n-norm', 'B-norm'
    
    line_n, = ax.plot(x, n, 'b')
    ax.set_xlabel(label_x)
    ax.set_ylabel(label_n, color='b')
    ax.tick_params(axis='y', labelcolor='b')
    ax.grid(True)

    ax2 = ax.twinx()
    line_b, = ax2.plot(x, B, 'r')
    ax2.set_ylabel(label_b, color='r')
    ax2.tick_params(axis='y', labelcolor='r')
    ax2.ticklabel_format(
        axis = 'y', 
        style='sci',scilimits = (-2,2)
    )  # !!! has to be after the subplots !!!
    ax.set_title(title)

    lines_dict = {
        name_line_n: line_n,
        name_line_B: line_b
    }
    return lines_dict

def find_cutoff_UHres(n, B, w_ref):
    nx_work = int(len(n)/2)
    n_work = n[0:nx_work]
    B_work = B[0:nx_work]

    wln = wl_func(n_work, B_work)/w_ref
    abs_diff = np.abs(wln-1)
    id_x_L_cutoff = np.where(abs_diff == np.min(abs_diff))[0][0] 

    wuhn = wuh_func(n_work, B_work)/w_ref
    abs_diff = np.abs(wuhn-1)
    id_x_UH_res = np.where(abs_diff == np.min(abs_diff))[0][0]

    return id_x_L_cutoff, id_x_UH_res

def plot_ws(ax, x, n, B, w_ref, wn):
    # x - norm. space grid (to r0)
    # n, B in Gauss units
    # w_ref - reference frequency in rad/s

    # find resonance and cutoff position:
    id_x_L_cutoff, id_x_UH_res = find_cutoff_UHres(n, B, w_ref)

    # line_wc,  = ax.plot(x, wc_func(B),     'black', linestyle='-', label='$\omega_c$')
    # line_wp,  = ax.plot(x, wp_func(n),     'grey',  linestyle='-', label='$\omega_p$')
    # line_wl,  = ax.plot(x, wl_func(n, B),  'b',     linestyle='-', label='$\omega_{L}$')
    # line_wuh, = ax.plot(x, wuh_func(n, B), 'g',     linestyle='-', label='$\omega_{uh}$')
    # line_wr,  = ax.plot(x, wr_func(n, B),  'r',     linestyle='-', label='$\omega_{R}$')
    # lw1 = ax.axhline(w_ref, color='b', linestyle = "--", linewidth = 2, label = "$k*c$")
    # lw_L_cutoff = ax.axvline(x[id_x_L_cutoff], color='b', linestyle = ":", linewidth = 2, label = "L-cutoff")
    # lw_UH_res = ax.axvline(x[id_x_UH_res], color='g', linestyle = ":", linewidth = 2, label = "UH-res")
    # ax.set_xlabel('$x/r_0$')
    # ax.set_ylabel('$\omega[rad/s]$')
    # ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left')
    # ax.grid(True)
    # plt.subplots_adjust(right=0.84)
    # ax.set_title("$Frequency\ profiles$")

    line_wc,  = ax.plot(x, wc_func(B)/wn,     'black', linestyle='-', label='$\omega_c$')
    line_wp,  = ax.plot(x, wp_func(n)/wn,     'grey',  linestyle='-', label='$\omega_p$')
    line_wl,  = ax.plot(x, wl_func(n, B)/wn,  'b',     linestyle='-', label='$\omega_{L}$')
    line_wuh, = ax.plot(x, wuh_func(n, B)/wn, 'g',     linestyle='-', label='$\omega_{uh}$')
    line_wr,  = ax.plot(x, wr_func(n, B)/wn,  'r',     linestyle='-', label='$\omega_{R}$')
    lw1 = ax.axhline(w_ref/wn, color='b', linestyle = "--", linewidth = 2, label = "$k*c$")
    lw_L_cutoff = ax.axvline(x[id_x_L_cutoff], color='b', linestyle = ":", linewidth = 2, label = "L-cutoff")
    lw_UH_res = ax.axvline(x[id_x_UH_res], color='g', linestyle = ":", linewidth = 2, label = "UH-res")
    ax.set_xlabel('$x/r_0$')
    ax.set_ylabel('$\omega/\omega_p$')
    ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left')
    ax.grid(True)
    plt.subplots_adjust(right=0.84)
    ax.set_title("$Frequency\ profiles$")

    lines_dict = {
        'w-ref': lw1, 'wc': line_wc, 'wp': line_wp, 'wl': line_wl, 'wuh': line_wuh, 'wr': line_wr,
        'L-cutoff': lw_L_cutoff, 'UH-res': lw_UH_res
    }
    return lines_dict, id_x_L_cutoff, id_x_UH_res

def w_grids(n, B):
    wc_grid = wc_func(B)
    wp_grid = wp_func(n)
    
    wl_grid  = wl_func(n, B)
    wuh_grid = wuh_func(n, B)
    wr_grid  = wr_func(n, B)

    return wc_grid, wp_grid, wl_grid, wuh_grid, wr_grid

def cma_grids(n, B, w_ref):
    _, _, wl_grid, wuh_grid, wr_grid = w_grids(n, B)
    
    wc2_ref_grid = wc_func(B)**2 / w_ref**2 
    wp2_ref_grid = wp_func(n)**2 / w_ref**2 

    wc2_l_grid = wc_func(B)**2 / wl_grid**2 
    wp2_l_grid = wp_func(n)**2 / wl_grid**2 

    wc2_r_grid = wc_func(B)**2 / wr_grid**2 
    wp2_r_grid = wp_func(n)**2 / wr_grid**2 

    wc2_uh_grid = wc_func(B)**2 / wuh_grid**2 
    wp2_uh_grid = wp_func(n)**2 / wuh_grid**2

    res = {
        'x-ref': wp2_ref_grid, 'y-ref': wc2_ref_grid,
        'x-l':   wp2_l_grid,   'y-l':   wc2_l_grid,
        'x-uh':  wp2_uh_grid,  'y-uh':  wc2_uh_grid,
        'x-r':   wp2_r_grid,   'y-r':   wc2_r_grid
    }
    return res

def plot_cma(ax, n, B, w_ref):
    # n, B in Gauss units
    # w_ref - reference frequency in rad/s
    gs = cma_grids(n, B, w_ref)

    line_cma_ref, = ax.plot(gs['x-ref'], gs['y-ref'], 'black', linestyle=':', label='$\omega_{ant}$')
    line_cma_l,   = ax.plot(gs['x-l'],   gs['y-l'],   'b',     linestyle='-', label='$\omega_{L}$')
    line_cma_uh,  = ax.plot(gs['x-uh'],  gs['y-uh'],  'g',     linestyle='-', label='$\omega_{uh}$')
    line_cma_r,   = ax.plot(gs['x-r'],   gs['y-r'],   'r',     linestyle='-', label='$\omega_{R}$')
    ax.set_xlabel('$\omega_p^2/\omega^2$')
    ax.set_ylabel('$\omega_c^2/\omega^2$')
    ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left')
    ax.grid(True)
    ax.set_title("$CMA\ diagram$")
    
    ax.set_xlim([0,2])
    ax.set_ylim([0,1])

    lines_dict = {
        'cma-ref': line_cma_ref, 'cma-l': line_cma_l, 'cma-uh': line_cma_uh, 'cma-r': line_cma_r
    }
    
    return lines_dict

def plot_tau_k(ax, x_plot, n, B, kx, dx, coef_t_norm, id_x_L_cutoff, id_x_UH_res):
    tau = kx * dx * len(x_plot) / wuh_func(n, B)
    tau_norm = tau / coef_t_norm

    line_tau, = ax.plot(x_plot, tau_norm, 'blue', linestyle='-', label='tau-k')
    lw_L_cutoff = ax.axvline(x_plot[id_x_L_cutoff], color='b', linestyle = ":", linewidth = 2, label = "L-cutoff")
    lw_UH_res = ax.axvline(x_plot[id_x_UH_res], color='g', linestyle = ":", linewidth = 2, label = "UH-res")

    ax.set_xlabel('$x/r_0$')
    ax.set_ylabel('$\\tau_k$')
    ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left')
    ax.grid(True)
    ax.set_title("$Time\ of\ k\ change$")

    lines_dict = {
        'tau-k': line_tau,
        'tauk-L-cutoff': lw_L_cutoff, 'tauk-UH-res': lw_UH_res
    }

    return lines_dict

def interact_nB(lines, x, x_plot, sel_norm, dx, kx, n0_input, x_n0_input, x_nw_input, B0_input, R0_input, kx_input):
    n = profile_n(n0_input, x, x_n0_input, x_nw_input)
    B = profile_B(B0_input, R0_input, x)

    w_ref_input = kx_input*c_light
    wc_grid, wp_grid, wl_grid, wuh_grid, wr_grid = w_grids(n, B)
    
    coef_v_norm, coef_t_norm, coef_r_norm, coef_e_norm, coef_n_norm = normalization(n0_input, B0_input, sel_norm) 
    n_norm = n / coef_n_norm
    B_norm = B / coef_e_norm

    tau = kx * dx * len(x) / wuh_func(n, B)
    tau_norm = tau / coef_t_norm

    wln = wl_func(n, B)/w_ref_input
    abs_diff = np.abs(wln-1)
    id_x_L_cutoff = np.where(abs_diff == np.min(abs_diff))[0][0]
    x_plot_L_cutoff = x_plot[id_x_L_cutoff] 

    wuhn = wuh_func(n, B)/w_ref_input
    abs_diff = np.abs(wuhn-1)
    id_x_UH_res = np.where(abs_diff == np.min(abs_diff))[0][0]
    x_plot_UH_res = x_plot[id_x_UH_res]

    gs = cma_grids(n, B, w_ref_input)
    
    print("n0 = {:0.3e} cm^(-3)".format(n0_input))
    print("B0 = {:0.3e} gauss".format(B0_input))
    print("R0 = {:0.1f} cm".format(R0_input))
    print('kx[cm^(-1)] = {:0.3f}'.format(kx_input))
    print('k*c[rad/s] = {:0.3e}'.format(w_ref_input))
    print('x/r0 of L-cutoff = {:0.3f}'.format(x_plot_L_cutoff))
    print('x/r0 of UH-res = {:0.3f}'.format(x_plot_UH_res))
    print('x/r0 diff between L-cutoff and UH-res = {:0.6f}'
        .format(np.abs(x_plot_L_cutoff - x_plot_UH_res)))
    
    lines['n'].set_ydata(n)
    lines['B'].set_ydata(B)
    
    lines['n-norm'].set_ydata(n_norm)
    lines['B-norm'].set_ydata(B_norm)
    
    lines['w-ref'].set_ydata(w_ref_input)
    lines['wc'].set_ydata(wc_grid)
    lines['wp'].set_ydata(wp_grid)
    lines['wl'].set_ydata(wl_grid)
    lines['wuh'].set_ydata(wuh_grid)
    lines['wr'].set_ydata(wr_grid)
    lines['L-cutoff'].set_xdata(x_plot_L_cutoff)
    lines['UH-res'].set_xdata(x_plot_UH_res)
    
    lines['cma-ref'].set_xdata(gs['x-ref'])
    lines['cma-ref'].set_ydata(gs['y-ref'])
    lines['cma-l'].set_xdata(gs['x-l'])
    lines['cma-l'].set_ydata(gs['y-l'])
    lines['cma-uh'].set_xdata(gs['x-uh']) 
    lines['cma-uh'].set_ydata(gs['y-uh'])
    lines['cma-r'].set_xdata(gs['x-r']) 
    lines['cma-r'].set_ydata(gs['y-r'])

    lines['tau-k'].set_ydata(tau_norm)
    lines['tauk-L-cutoff'].set_xdata(x_plot_L_cutoff)
    lines['tauk-UH-res'].set_xdata(x_plot_UH_res)
    return 

def plot_nB_w_cma_profiles(dd):
    x = dd['x']
    x_plot = dd['x-plot']
    n = dd['n']
    B = dd['B']
    w_ref = dd['w-ref']
    w_norm = 1./ dd["coef-t-norm"]
    kx = dd['kx']
    dx = dd['dx']
    coef_t_norm = dd['coef-t-norm']
    sel_norm = dd['sel-norm']

    lines_dict = {}
    fig, axs = plt.subplots(3,2, figsize=(9.6,9))

    lines_new = plot_nB(axs[0,0], x, n, B, False) 
    lines_dict.update(lines_new)

    lines_new = plot_nB(axs[0,1], x_plot, dd['n-norm'], dd['B-norm'], True)
    lines_dict.update(lines_new)

    lines_new, id_x_L_cutoff, id_x_UH_res = plot_ws(axs[1,0], x_plot, n, B, w_ref, w_norm)
    lines_dict.update(lines_new)

    lines_new = plot_cma(axs[1,1], n, B, w_ref)
    lines_dict.update(lines_new)

    line_tau = plot_tau_k(axs[2,0], x_plot, n, B, kx, dx, coef_t_norm, id_x_L_cutoff, id_x_UH_res)
    lines_dict.update(line_tau)

    plt.tight_layout()
    plt.show()

    dd_int = dd['interact']
    n_int = dd_int['n']
    x_n0_int = dd_int['x-n0']
    x_nw_int = dd_int['x-nw']
    B_int = dd_int['B']
    R_int = dd_int['R']
    kx_int = dd_int['kx']

    xx_n = interact(
        lambda n0_input, x_n0_input, x_nw_input, B0_input, R0_input, kx_input: \
            interact_nB(lines_dict, x, x_plot, sel_norm, dx, kx,
                n0_input, x_n0_input, x_nw_input, B0_input, R0_input, kx_input),
        
        n0_input=widgets.FloatSlider(
            min=n_int[1], 
            max=n_int[2], 
            value=n_int[0], 
            step=n_int[3],
            description='n0 [cm^(-3)]:',
            readout_format='.3e'
        ),
        x_n0_input=widgets.FloatSlider(
            min=x_n0_int[1], 
            max=x_n0_int[2], 
            value=x_n0_int[0], 
            step=x_n0_int[3],
            description='x of n0 [cm]:',
            readout_format='.3e'
        ),
        x_nw_input=widgets.FloatSlider(
            min=x_nw_int[1], 
            max=x_nw_int[2], 
            value=x_nw_int[0], 
            step=x_nw_int[3],
            description='width of n [cm]:',
            readout_format='.3e'
        ),
        B0_input=widgets.FloatSlider(
            min=B_int[1], 
            max=B_int[2], 
            value=B_int[0], 
            step=B_int[3],
            description='B0 [gauss]:',
            readout_format='.3e'
        ),
        R0_input=widgets.FloatSlider(
            min=R_int[1], 
            max=R_int[2], 
            value=R_int[0], 
            step=R_int[3],
            description='R0 [cm]:',
            readout_format='.1f'
        ),
        kx_input=widgets.FloatSlider(
            min=kx_int[1], 
            max=kx_int[2], 
            value=kx_int[0], 
            step=kx_int[3],
            description='kx [cm^(-1)]:',
            readout_format='.3f'
        ),
    )

    res = {
        'id-x-L-cutoff': id_x_L_cutoff,
        'id-x-UH-res': id_x_UH_res
    }
    return res


def interact_nB_outgoing(
    lines, x, x_plot, sel_norm, dx, kx, 
    n0_input, x_n0_input, x_nw_input, B0_input, R0_input, 
    n0_input_2, x_n0_input_2, x_nw_input_2, B0_input_2, R0_input_2, x_B_left_input,
    kx_input, n_gap
):
    n = profile_n(n0_input, x, x_n0_input, x_nw_input, n0_input_2, x_n0_input_2, x_nw_input_2)
    B = profile_B(B0_input, R0_input, x, B0_input_2, R0_input_2, x_B_left_input, n_gap)

    w_ref_input = kx_input*c_light
    wc_grid, wp_grid, wl_grid, wuh_grid, wr_grid = w_grids(n, B)
    
    coef_v_norm, coef_t_norm, coef_r_norm, coef_e_norm, coef_n_norm = normalization(n0_input, B0_input, sel_norm) 
    n_norm = n / coef_n_norm
    B_norm = B / coef_e_norm

    w_norm = 1./coef_t_norm

    tau = kx * dx * len(x) / wuh_func(n, B)
    tau_norm = tau / coef_t_norm

    id_x_L_cutoff, id_x_UH_res = find_cutoff_UHres(n, B, w_ref_input)
    x_plot_L_cutoff = x_plot[id_x_L_cutoff]
    x_plot_UH_res   = x_plot[id_x_UH_res] 

    # wln = wl_func(n, B)/w_ref_input
    # abs_diff = np.abs(wln-1)
    # id_x_L_cutoff = np.where(abs_diff == np.min(abs_diff))[0][0]
    # x_plot_L_cutoff = x_plot[id_x_L_cutoff] 

    # wuhn = wuh_func(n, B)/w_ref_input
    # abs_diff = np.abs(wuhn-1)
    # id_x_UH_res = np.where(abs_diff == np.min(abs_diff))[0][0]
    # x_plot_UH_res = x_plot[id_x_UH_res]

    gs = cma_grids(n, B, w_ref_input)
    
    print("n0 = {:0.3e} cm^(-3)".format(n0_input))
    print("B0 = {:0.3e} gauss".format(B0_input))
    print("R0 = {:0.1f} cm".format(R0_input))
    print('kx[cm^(-1)] = {:0.3f}'.format(kx_input))
    print('k*c[rad/s] = {:0.3e}'.format(w_ref_input))
    print('x/r0 of L-cutoff = {:0.3f}'.format(x_plot_L_cutoff))
    print('x/r0 of UH-res = {:0.3f}'.format(x_plot_UH_res))
    print('x/r0 diff between L-cutoff and UH-res = {:0.6f}'
        .format(np.abs(x_plot_L_cutoff - x_plot_UH_res)))
    
    lines['n'].set_ydata(n)
    lines['B'].set_ydata(B)
    
    lines['n-norm'].set_ydata(n_norm)
    lines['B-norm'].set_ydata(B_norm)
    
    # lines['w-ref'].set_ydata(w_ref_input)
    # lines['wc'].set_ydata(wc_grid)
    # lines['wp'].set_ydata(wp_grid)
    # lines['wl'].set_ydata(wl_grid)
    # lines['wuh'].set_ydata(wuh_grid)
    # lines['wr'].set_ydata(wr_grid)
    lines['w-ref'].set_ydata(w_ref_input/w_norm)
    lines['wc'].set_ydata(wc_grid/w_norm)
    lines['wp'].set_ydata(wp_grid/w_norm)
    lines['wl'].set_ydata(wl_grid/w_norm)
    lines['wuh'].set_ydata(wuh_grid/w_norm)
    lines['wr'].set_ydata(wr_grid/w_norm)
    lines['L-cutoff'].set_xdata(x_plot_L_cutoff)
    lines['UH-res'].set_xdata(x_plot_UH_res)
    
    lines['cma-ref'].set_xdata(gs['x-ref'])
    lines['cma-ref'].set_ydata(gs['y-ref'])
    lines['cma-l'].set_xdata(gs['x-l'])
    lines['cma-l'].set_ydata(gs['y-l'])
    lines['cma-uh'].set_xdata(gs['x-uh']) 
    lines['cma-uh'].set_ydata(gs['y-uh'])
    lines['cma-r'].set_xdata(gs['x-r']) 
    lines['cma-r'].set_ydata(gs['y-r'])

    lines['tau-k'].set_ydata(tau_norm)
    lines['tauk-L-cutoff'].set_xdata(x_plot_L_cutoff)
    lines['tauk-UH-res'].set_xdata(x_plot_UH_res)
    return 

def plot_nB_w_cma_profiles_outgoing(dd):
    x = dd['x']
    x_plot = dd['x-plot']
    n = dd['n']
    B = dd['B']
    w_ref = dd['w-ref']
    w_norm = 1./ dd["coef-t-norm"]
    kx = dd['kx']
    dx = dd['dx']
    coef_t_norm = dd['coef-t-norm']
    sel_norm = dd['sel-norm']

    lines_dict = {}
    fig, axs = plt.subplots(3,2, figsize=(9.6,9))

    lines_new = plot_nB(axs[0,0], x, n, B, False) 
    lines_dict.update(lines_new)

    lines_new = plot_nB(axs[0,1], x_plot, dd['n-norm'], dd['B-norm'], True)
    lines_dict.update(lines_new)

    lines_new, id_x_L_cutoff, id_x_UH_res = plot_ws(axs[1,0], x_plot, n, B, w_ref, w_norm)
    lines_dict.update(lines_new)

    lines_new = plot_cma(axs[1,1], n, B, w_ref)
    lines_dict.update(lines_new)

    line_tau = plot_tau_k(axs[2,0], x_plot, n, B, kx, dx, coef_t_norm, id_x_L_cutoff, id_x_UH_res)
    lines_dict.update(line_tau)

    plt.tight_layout()
    plt.show()

    dd_int = dd['interact']
    n_int = dd_int['n']
    x_n0_int = dd_int['x-n0']
    x_nw_int = dd_int['x-nw']
    B_int = dd_int['B']
    R_int = dd_int['R']

    n_int_2 = dd_int['n-2']
    x_n0_int_2 = dd_int['x-n0-2']
    x_nw_int_2 = dd_int['x-nw-2']
    B_int_2 = dd_int['B-2']
    R_int_2 = dd_int['R-2']
    x_B_left_int = dd_int['x-B-left']

    kx_int = dd_int['kx']
    n_gap = dd["np-gap-B"]

    xx_n = interact(
        lambda n0_input, x_n0_input, x_nw_input, B0_input, R0_input, n0_input_2, x_n0_input_2, x_nw_input_2, B0_input_2, R0_input_2, x_B_left_input, kx_input: \
            interact_nB_outgoing(lines_dict, x, x_plot, sel_norm, dx, kx,
                n0_input,   x_n0_input,   x_nw_input,   B0_input,   R0_input, 
                n0_input_2, x_n0_input_2, x_nw_input_2, B0_input_2, R0_input_2, x_B_left_input,
                kx_input, n_gap),
        
        n0_input=widgets.FloatSlider(
            min=n_int[1], 
            max=n_int[2], 
            value=n_int[0], 
            step=n_int[3],
            description='n0 [cm^(-3)]:',
            readout_format='.3e'
        ),
        x_n0_input=widgets.FloatSlider(
            min=x_n0_int[1], 
            max=x_n0_int[2], 
            value=x_n0_int[0], 
            step=x_n0_int[3],
            description='x of n0 [cm]:',
            readout_format='.3e'
        ),
        x_nw_input=widgets.FloatSlider(
            min=x_nw_int[1], 
            max=x_nw_int[2], 
            value=x_nw_int[0], 
            step=x_nw_int[3],
            description='width of n [cm]:',
            readout_format='.3e'
        ),
        B0_input=widgets.FloatSlider(
            min=B_int[1], 
            max=B_int[2], 
            value=B_int[0], 
            step=B_int[3],
            description='B0 [gauss]:',
            readout_format='.3e'
        ),
        R0_input=widgets.FloatSlider(
            min=R_int[1], 
            max=R_int[2], 
            value=R_int[0], 
            step=R_int[3],
            description='R0 [cm]:',
            readout_format='.1f'
        ),

        n0_input_2=widgets.FloatSlider(
            min=n_int_2[1], 
            max=n_int_2[2], 
            value=n_int_2[0], 
            step=n_int_2[3],
            description='n0-2 [cm^(-3)]:',
            readout_format='.3e'
        ),
        x_n0_input_2=widgets.FloatSlider(
            min=x_n0_int_2[1], 
            max=x_n0_int_2[2], 
            value=x_n0_int_2[0], 
            step=x_n0_int_2[3],
            description='x of n0-2 [cm]:',
            readout_format='.3e'
        ),
        x_nw_input_2=widgets.FloatSlider(
            min=x_nw_int_2[1], 
            max=x_nw_int_2[2], 
            value=x_nw_int_2[0], 
            step=x_nw_int_2[3],
            description='width of n-2 [cm]:',
            readout_format='.3e'
        ),
        B0_input_2=widgets.FloatSlider(
            min=B_int_2[1], 
            max=B_int_2[2], 
            value=B_int_2[0], 
            step=B_int_2[3],
            description='B0-2 [gauss]:',
            readout_format='.3e'
        ),
        R0_input_2=widgets.FloatSlider(
            min=R_int_2[1], 
            max=R_int_2[2], 
            value=R_int_2[0], 
            step=R_int_2[3],
            description='R0-2 [cm]:',
            readout_format='.1f'
        ),
        x_B_left_input=widgets.FloatSlider(
            min=x_B_left_int[1], 
            max=x_B_left_int[2], 
            value=x_B_left_int[0], 
            step=x_B_left_int[3],
            description='x-B-left [cm]:',
            readout_format='.1f'
        ),

        kx_input=widgets.FloatSlider(
            min=kx_int[1], 
            max=kx_int[2], 
            value=kx_int[0], 
            step=kx_int[3],
            description='kx [cm^(-1)]:',
            readout_format='.3f'
        ),
    )

    res = {
        'id-x-L-cutoff': id_x_L_cutoff,
        'id-x-UH-res': id_x_UH_res
    }
    return res


def interact_profiles(
    lines, x, x_plot, sel_norm, dx, kx, 
    n0_input, x_n0_input, x_nw_input, B0_input, R0_input, 
    kx_input
):
    n = profile_n(n0_input, x, x_n0_input, x_nw_input)
    B = profile_B(B0_input, R0_input, x)

    w_ref_input = kx_input*c_light
    wc_grid, wp_grid, wl_grid, wuh_grid, wr_grid = w_grids(n, B)
    
    coef_v_norm, coef_t_norm, coef_r_norm, coef_e_norm, coef_n_norm = normalization(n0_input, B0_input, sel_norm) 
    n_norm = n / coef_n_norm
    B_norm = B / coef_e_norm

    tau = kx * dx * len(x) / wuh_func(n, B)
    tau_norm = tau / coef_t_norm

    id_x_L_cutoff, id_x_UH_res = find_cutoff_UHres(n, B, w_ref_input)
    x_plot_L_cutoff = x_plot[id_x_L_cutoff]
    x_plot_UH_res   = x_plot[id_x_UH_res] 

    # wln = wl_func(n, B)/w_ref_input
    # abs_diff = np.abs(wln-1)
    # id_x_L_cutoff = np.where(abs_diff == np.min(abs_diff))[0][0]
    # x_plot_L_cutoff = x_plot[id_x_L_cutoff] 

    # wuhn = wuh_func(n, B)/w_ref_input
    # abs_diff = np.abs(wuhn-1)
    # id_x_UH_res = np.where(abs_diff == np.min(abs_diff))[0][0]
    # x_plot_UH_res = x_plot[id_x_UH_res]

    gs = cma_grids(n, B, w_ref_input)
    
    print("n0 = {:0.3e} cm^(-3)".format(n0_input))
    print("B0 = {:0.3e} gauss".format(B0_input))
    print("R0 = {:0.1f} cm".format(R0_input))
    print('kx[cm^(-1)] = {:0.3f}'.format(kx_input))
    print('k*c[rad/s] = {:0.3e}'.format(w_ref_input))
    print('x/r0 of L-cutoff = {:0.3f}'.format(x_plot_L_cutoff))
    print('x/r0 of UH-res = {:0.3f}'.format(x_plot_UH_res))
    print('x/r0 diff between L-cutoff and UH-res = {:0.6f}'
        .format(np.abs(x_plot_L_cutoff - x_plot_UH_res)))
    
    lines['n'].set_ydata(n)
    lines['B'].set_ydata(B)
    
    lines['n-norm'].set_ydata(n_norm)
    lines['B-norm'].set_ydata(B_norm)
    
    lines['w-ref'].set_ydata(w_ref_input)
    lines['wc'].set_ydata(wc_grid)
    lines['wp'].set_ydata(wp_grid)
    lines['wl'].set_ydata(wl_grid)
    lines['wuh'].set_ydata(wuh_grid)
    lines['wr'].set_ydata(wr_grid)
    lines['L-cutoff'].set_xdata(x_plot_L_cutoff)
    lines['UH-res'].set_xdata(x_plot_UH_res)
    
    lines['cma-ref'].set_xdata(gs['x-ref'])
    lines['cma-ref'].set_ydata(gs['y-ref'])
    lines['cma-l'].set_xdata(gs['x-l'])
    lines['cma-l'].set_ydata(gs['y-l'])
    lines['cma-uh'].set_xdata(gs['x-uh']) 
    lines['cma-uh'].set_ydata(gs['y-uh'])
    lines['cma-r'].set_xdata(gs['x-r']) 
    lines['cma-r'].set_ydata(gs['y-r'])

    lines['tau-k'].set_ydata(tau_norm)
    lines['tauk-L-cutoff'].set_xdata(x_plot_L_cutoff)
    lines['tauk-UH-res'].set_xdata(x_plot_UH_res)
    return 


def plot_profiles(dd):
    x = dd['x']
    x_plot = dd['x-plot']
    n = dd['n']
    B = dd['B']
    w_ref = dd['w-ref']
    w_norm = 1./ dd["coef-t-norm"]
    kx = dd['kx']
    dx = dd['dx']
    coef_t_norm = dd['coef-t-norm']
    sel_norm = dd['sel-norm']

    lines_dict = {}
    fig, axs = plt.subplots(3,2, figsize=(9.6,9))

    lines_new = plot_nB(axs[0,0], x, n, B, False) 
    lines_dict.update(lines_new)

    lines_new = plot_nB(axs[0,1], x_plot, dd['n-norm'], dd['B-norm'], True)
    lines_dict.update(lines_new)

    lines_new, id_x_L_cutoff, id_x_UH_res = plot_ws(axs[1,0], x_plot, n, B, w_ref, w_norm)
    lines_dict.update(lines_new)

    lines_new = plot_cma(axs[1,1], n, B, w_ref)
    lines_dict.update(lines_new)

    line_tau = plot_tau_k(axs[2,0], x_plot, n, B, kx, dx, coef_t_norm, id_x_L_cutoff, id_x_UH_res)
    lines_dict.update(line_tau)

    plt.tight_layout()
    plt.show()

    dd_int = dd['interact']
    n_int = dd_int['n']
    x_n0_int = dd_int['x-n0']
    x_nw_int = dd_int['x-nw']
    B_int = dd_int['B']
    R_int = dd_int['R']

    kx_int = dd_int['kx']

    xx_n = interact(
        lambda n0_input, x_n0_input, x_nw_input, B0_input, R0_input, kx_input: \
            interact_profiles(lines_dict, x, x_plot, sel_norm, dx, kx,
                n0_input,   x_n0_input,   x_nw_input,   B0_input,   R0_input, 
                kx_input),
        
        n0_input=widgets.FloatSlider(
            min=n_int[1], 
            max=n_int[2], 
            value=n_int[0], 
            step=n_int[3],
            description='n0 [cm^(-3)]:',
            readout_format='.3e'
        ),
        x_n0_input=widgets.FloatSlider(
            min=x_n0_int[1], 
            max=x_n0_int[2], 
            value=x_n0_int[0], 
            step=x_n0_int[3],
            description='x of n0 [cm]:',
            readout_format='.3e'
        ),
        x_nw_input=widgets.FloatSlider(
            min=x_nw_int[1], 
            max=x_nw_int[2], 
            value=x_nw_int[0], 
            step=x_nw_int[3],
            description='width of n [cm]:',
            readout_format='.3e'
        ),
        B0_input=widgets.FloatSlider(
            min=B_int[1], 
            max=B_int[2], 
            value=B_int[0], 
            step=B_int[3],
            description='B0 [gauss]:',
            readout_format='.3e'
        ),
        R0_input=widgets.FloatSlider(
            min=R_int[1], 
            max=R_int[2], 
            value=R_int[0], 
            step=R_int[3],
            description='R0 [cm]:',
            readout_format='.1f'
        ),
        kx_input=widgets.FloatSlider(
            min=kx_int[1], 
            max=kx_int[2], 
            value=kx_int[0], 
            step=kx_int[3],
            description='kx [cm^(-1)]:',
            readout_format='.3f'
        ),
    )

    res = {
        'id-x-L-cutoff': id_x_L_cutoff,
        'id-x-UH-res': id_x_UH_res
    }
    return res


def set_init_variables_Q_outgoing(Q0, x0, dd):
    sel_norm = dd['sel-norm']
    if sel_norm != 'c':
        print('error: wrong normalization.')
        exit(-1)

    w_ref_norm = dd['w-ref-norm']
    nx, nt = dd['nx'], dd['nt']
    tau = dd['dt-norm']
    x_plot = dd['x-plot']

    Q = np.zeros((nx,nt), dtype = np.complex)
    vx = np.zeros((nx,nt), dtype = np.complex)
    vy = np.zeros((nx,nt), dtype = np.complex)
    Ex = np.zeros((nx,nt), dtype = np.complex)
    Ey = np.zeros((nx,nt), dtype = np.complex)
    Bz = np.zeros((nx,nt), dtype = np.complex)


    id_Q_1, _ = mix.find(x_plot, x0)
    id_Q_2 = id_Q_1 - 1

    Q[id_Q_1,0] = Q0 * np.exp(-1j * w_ref_norm * tau)
    Q[id_Q_1,1] = Q0

    Q[id_Q_2,0] = Q0 * np.exp(-1j * w_ref_norm * tau)
    Q[id_Q_2,1] = Q0

    res = {
        'vx': vx, 'vy': vy, 'Ex': Ex, 'Ey': Ey, 'Bz': Bz, 'Q': Q, 'ids-Q': [id_Q_1, id_Q_2]
    }
    return res

def compute_Q_outgoing_DD(coef_coupling, vars_data, dd):
    print("Computation...", end='\r')

    sel_norm = dd['sel-norm']
    if sel_norm != 'c':
        print('error: wrong normalization.')
        exit(-1)

    nx, nt = dd['nx'], dd['nt']
    tau, h = dd['dt-norm'], dd['dx-norm']
    t_norm = dd['t-norm']
    n_norm = dd['n-norm']
    Bb_norm = dd['B-norm']
    w_ref_norm = dd['w-ref-norm']
    
    vx = np.array(vars_data['vx'])
    vy = np.array(vars_data['vy'])
    Ex = np.array(vars_data['Ex'])
    Ey = np.array(vars_data['Ey'])
    Bz = np.array(vars_data['Bz'])
    Q = np.array(vars_data['Q'])

    ids_Q = vars_data['ids-Q']

    # REMARK: there are not ghost numbers in the x-grid
    def one_time_step(
        tau, h, t,  
        vx_prev, vy_prev, Ex_prev, Ey_prev, Bz_prev, Q_prev,
        vx, vy, Ex, Ey, Bz, Q
    ):
        tau2 = 2 * tau

        vx_new = np.zeros(nx, dtype = np.complex)
        vy_new = np.zeros(nx, dtype = np.complex)
        Ex_new = np.zeros(nx, dtype = np.complex)
        Ey_new = np.zeros(nx, dtype = np.complex)
        Bz_new = np.zeros(nx, dtype = np.complex)
        Q_new  = np.zeros(nx, dtype = np.complex)

        for id_Q in ids_Q:
            Q_new[id_Q] = Q_prev[id_Q] + tau2 * 1j*w_ref_norm * Q[id_Q] + 1j*tau2 * coef_coupling * Bz[id_Q]

        for i in range(nx):
            temp = tau2 * np.sqrt(n_norm[i]) * Ex[i]
            vx_new[i] = vx_prev[i] - tau2 * Bb_norm[i] * vy[i] - temp
            
            temp = tau2 * np.sqrt(n_norm[i]) * Ey[i]
            vy_new[i] = vy_prev[i] + tau2 * Bb_norm[i] * vx[i] - temp
                                            
            Ex_new[i] = Ex_prev[i] + tau2*np.sqrt(n_norm[i]) * vx[i]

            if i > 0 and i < (nx-1):
                temp = tau/h * (Bz[i+1] - Bz[i-1])
                Ey_new[i] = Ey_prev[i] + tau2*np.sqrt(n_norm[i]) * vy[i] - temp
                Bz_new[i] = Bz_prev[i] - tau/h * (Ey[i+1] - Ey[i-1]) + 1j*tau2 * coef_coupling * Q[i]
                
        return vx_new, vy_new, Ex_new, Ey_new, Bz_new, Q_new

    # *** computation ***
    for it in range(1, nt-1):
        if it%100 == 0:
            print("Computation... id-t = {:d}".format(it), end ='\r')

        vx[:,it+1], vy[:,it+1], Ex[:,it+1], Ey[:,it+1], Bz[:,it+1], Q[:,it+1] = \
            one_time_step(
                tau, h, t_norm[it], 
                vx[:,it-1], vy[:,it-1], Ex[:,it-1], Ey[:,it-1], Bz[:,it-1], Q[:,it-1],
                vx[:,it],   vy[:,it],   Ex[:,it],   Ey[:,it],   Bz[:,it], Q[:,it]
            )
    print()

    # energy density
    Wv = np.abs(vx)**2 + np.abs(vy)**2
    We = np.abs(Ex)**2 + np.abs(Ey)**2
    Wb = np.abs(Bz)**2
    Wq = np.abs(Q)**2
    W = Wv + We + Wb + Wq
        
    print("Done")
    print("(sumx W[1] - sumx W[-1])/sumx W[1]: {:0.3e}".format( 
        (np.sum(W[:,1]) - np.sum(W[:,-1])) / np.sum(W[:,1]) 
    ))

    res = {
        'vx': vx, 'vy': vy, 'Ex': Ex, 'Ey': Ey, 'Bz': Bz, 'Q': Q,
        'W': W, 'Wv': Wv, 'We': We, 'Wb': Wb, 'Wq': Wq
    }
    return res