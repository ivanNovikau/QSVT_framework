import numpy as np

import pylib.mix as mix
import pylib.ymath as ymath
import pylib.Global_variables as GLO
import pylib.plib as plib
import pylib.measurement as mse

import XW_waves.fluid_waves_1d as fw

def reload():
    mix.reload_module(mix)
    mix.reload_module(ymath)
    mix.reload_module(GLO)
    mix.reload_module(plib)
    mix.reload_module(mse)
    mix.reload_module(fw)



def init(case_to_model):
    if mix.compare_two_strings(case_to_model, "tiny-n3"):
        dd = model_tiny_case_n3()
    elif mix.compare_two_strings(case_to_model, "tiny-n4"):
        dd = model_tiny_case_n4()
    elif mix.compare_two_strings(case_to_model, "tiny-n5"):
        dd = model_tiny_case_n5()
    elif mix.compare_two_strings(case_to_model, "small"):
        dd = model_small_case()
    elif mix.compare_two_strings(case_to_model, "ref"):
        dd = model_ref_case()
    else:
        print("There is not a case with a name ", case_to_model)
        return

    # --- initial source value ---
    dd["Q0-qsp"] = 1. / np.sqrt(2)
    beta = 1./ dd["norm-of-h"] 
    dd["Q0"] = dd["Q0-qsp"] / beta

    # --- print Classical parameters ---
    print("---------------------------------------------------------------------")
    print("--- Classical parameters ---")
    print("nqx, nx = {:d}, {:d}".format(dd["nqx"], dd["nx"]))
    print("t, nt = {:0.3f}, {:d}".format(dd["t0"], dd["nt"]))
    print("norm. time step: {:0.3e}".format(dd["dt-norm"]))
    print("Courant number is: {:0.3e}".format(dd["Courant"]))

    # --- print QC data ---
    print("\n---------------------------------------------------------------------")
    print("--- QC parameters ---")
    print("QC time step: {:0.3f}".format(dd["dt-norm"] * dd["norm-of-h"]))
    print("QC normalization of a Hamiltonian matrix: {:0.3f}".format(dd["norm-of-h"]))
    print("QC 1/norm_of_H: {:0.3e}".format(1./dd["norm-of-h"]))
    print("QC source init. value: {:0.3e}".format(dd["Q0-qsp"]))

    qsp_t_desired = dd["t0"] * dd["norm-of-h"] 
    nt_resulting_qsp = int(dd["t0"] / dd["dt-norm"])
    print("QC time to simulate: {:0.3f}".format(qsp_t_desired))
    print("QC number of time steps: {:d}".format(nt_resulting_qsp))

    return dd


def model_tiny_case_n3():
    print("Simulation of a --- TINY --- case\n")

    path_root_classical = "./QSP/examples/XW_1D/classic-results/"
    path_qc = "./QSP/examples/XW_1D/tiny-case-n3/"
    fname_start = "small-"

    nqx = 3

    r0, nx = 2,  2**nqx   
    t0, nt = 100, 400
    n0, x_n0, x_nw = 2e13, -0.99 * r0, 0.20 * r0
    B0, R0 = 1e3, 167
    n0_2, x_n0_2, x_nw_2 = 1e12, 0.90 * r0, 0.18 * r0
    B0_2, R0_2, x_B_left = 6e3, 10, 0.4 * r0
    kx = 3.228
    x0 = -0.1  # where Q sits
    coef_coupling = 0.100 # coupling between Q and Bz

    np_gap_B = 3

    # coefficient to take into account additional multipliers in an oracle due to superpositions:
    # in 1-D fluid plasma X-wave problem: there is a superposition of states -> max multiplier is 2
    coef_superposition = 2 # might change for a different discretization scheme

    # --- prepare profiles and normalization ---
    dd_init = {
        "path-classical": path_root_classical + "/" + fname_start,
        "path-qc": path_qc,
        "r0": r0, "t0": t0,
        "n0": n0, "B0": B0, "R0": R0, "kx": kx,
        "nqx": nqx, "nx": nx, "nt": nt,
        "x_n0": x_n0, "x_nw": x_nw,
        "n0_2": n0_2, "x_n0_2": x_n0_2, "x_nw_2": x_nw_2,
        "B0_2": B0_2, "R0_2": R0_2, "x_B_left": x_B_left,
        "sel-norm": "c",
        "coef-coupling": coef_coupling,
        "coef-superposition": coef_superposition,
        "x-source": x0,
        "np-gap-B": np_gap_B
    }
    dd = fw.prepare_norm_profiles_outgoing(dd_init)
    return dd


def model_tiny_case_n4():
    print("Simulation of a --- TINY --- case\n")

    path_root_classical = "./QSP/examples/XW_1D/classic-results/"
    path_qc = "./QSP/examples/XW_1D/tiny-case-n4/"
    fname_start = "small-"

    nqx = 4

    r0, nx = 2,  2**nqx   
    t0, nt = 100, 400
    n0, x_n0, x_nw = 2e13, -0.99 * r0, 0.20 * r0
    B0, R0 = 1e3, 167
    n0_2, x_n0_2, x_nw_2 = 1e12, 0.90 * r0, 0.18 * r0
    B0_2, R0_2, x_B_left = 6e3, 10, 0.4 * r0
    kx = 3.228
    x0 = -0.1  # where Q sits
    coef_coupling = 0.100 # coupling between Q and Bz

    np_gap_B = 4

    # coefficient to take into account additional multipliers in an oracle due to superpositions:
    # in 1-D fluid plasma X-wave problem: there is a superposition of states -> max multiplier is 2
    coef_superposition = 2 # might change for a different discretization scheme

    # --- prepare profiles and normalization ---
    dd_init = {
        "path-classical": path_root_classical + "/" + fname_start,
        "path-qc": path_qc,
        "r0": r0, "t0": t0,
        "n0": n0, "B0": B0, "R0": R0, "kx": kx,
        "nqx": nqx, "nx": nx, "nt": nt,
        "x_n0": x_n0, "x_nw": x_nw,
        "n0_2": n0_2, "x_n0_2": x_n0_2, "x_nw_2": x_nw_2,
        "B0_2": B0_2, "R0_2": R0_2, "x_B_left": x_B_left,
        "sel-norm": "c",
        "coef-coupling": coef_coupling,
        "coef-superposition": coef_superposition,
        "x-source": x0,
        "np-gap-B": np_gap_B
    }
    dd = fw.prepare_norm_profiles_outgoing(dd_init)
    return dd


def model_tiny_case_n5():
    print("Simulation of a --- TINY --- case\n")

    path_root_classical = "./QSP/examples/XW_1D/classic-results/"
    path_qc = "./QSP/examples/XW_1D/tiny-case-n5/"
    fname_start = "small-"

    nqx = 5

    r0, nx = 2,  2**nqx   
    t0, nt = 100, 400
    n0, x_n0, x_nw = 2e13, -0.99 * r0, 0.20 * r0
    B0, R0 = 1e3, 167
    n0_2, x_n0_2, x_nw_2 = 1e12, 0.90 * r0, 0.18 * r0
    B0_2, R0_2, x_B_left = 6e3, 10, 0.4 * r0
    kx = 3.228
    x0 = -0.1  # where Q sits
    coef_coupling = 0.100 # coupling between Q and Bz

    np_gap_B = 5

    # coefficient to take into account additional multipliers in an oracle due to superpositions:
    # in 1-D fluid plasma X-wave problem: there is a superposition of states -> max multiplier is 2
    coef_superposition = 2 # might change for a different discretization scheme

    # --- prepare profiles and normalization ---
    dd_init = {
        "path-classical": path_root_classical + "/" + fname_start,
        "path-qc": path_qc,
        "r0": r0, "t0": t0,
        "n0": n0, "B0": B0, "R0": R0, "kx": kx,
        "nqx": nqx, "nx": nx, "nt": nt,
        "x_n0": x_n0, "x_nw": x_nw,
        "n0_2": n0_2, "x_n0_2": x_n0_2, "x_nw_2": x_nw_2,
        "B0_2": B0_2, "R0_2": R0_2, "x_B_left": x_B_left,
        "sel-norm": "c",
        "coef-coupling": coef_coupling,
        "coef-superposition": coef_superposition,
        "x-source": x0,
        "np-gap-B": np_gap_B
    }
    dd = fw.prepare_norm_profiles_outgoing(dd_init)
    return dd


def model_small_case():
    print("Simulation of a --- SMALL --- case\n")

    path_root_classical = "./QSP/examples/XW_1D/classic-results/"
    path_qc = "./QSP/examples/XW_1D/small-case/"
    fname_start = "small-"

    nqx = 6

    r0, nx = 2,  2**nqx   
    t0, nt = 100, 400
    n0, x_n0, x_nw = 2e13, -0.99 * r0, 0.20 * r0
    B0, R0 = 1e3, 167
    n0_2, x_n0_2, x_nw_2 = 1e12, 0.90 * r0, 0.18 * r0
    B0_2, R0_2, x_B_left = 6e3, 10, 0.4 * r0
    kx = 3.228
    x0 = -0.1  # where Q sits
    coef_coupling = 0.100 # coupling between Q and Bz

    np_gap_B = 10

    # coefficient to take into account additional multipliers in an oracle due to superpositions:
    # in 1-D fluid plasma X-wave problem: there is a superposition of states -> max multiplier is 2
    coef_superposition = 2 # might change for a different discretization scheme

    # --- prepare profiles and normalization ---
    dd_init = {
        "path-classical": path_root_classical + "/" + fname_start,
        "path-qc": path_qc,
        "r0": r0, "t0": t0,
        "n0": n0, "B0": B0, "R0": R0, "kx": kx,
        "nqx": nqx, "nx": nx, "nt": nt,
        "x_n0": x_n0, "x_nw": x_nw,
        "n0_2": n0_2, "x_n0_2": x_n0_2, "x_nw_2": x_nw_2,
        "B0_2": B0_2, "R0_2": R0_2, "x_B_left": x_B_left,
        "sel-norm": "c",
        "coef-coupling": coef_coupling,
        "coef-superposition": coef_superposition,
        "x-source": x0,
        "np-gap-B": np_gap_B
    }
    dd = fw.prepare_norm_profiles_outgoing(dd_init)
    return dd


def model_ref_case():
    print("Simulation of a --- REF --- case\n")

    path_root_classical = "./QSP/examples/XW_1D/classic-results/"
    path_qc = "./QSP/examples/XW_1D/ref-case/"
    fname_start = "ref-"

    nqx = 10

    r0, nx = 20,  2**nqx   
    t0, nt = 325, 1300
    n0, x_n0, x_nw = 2e13, -0.99 * r0, 0.20 * r0
    B0, R0 = 1e3, 167
    n0_2, x_n0_2, x_nw_2 = 1e12, 0.90 * r0, 0.18 * r0
    B0_2, R0_2, x_B_left = 6e3, 10, 0.4 * r0
    kx = 3.228
    x0 = -0.1  # where Q sits
    coef_coupling = 0.100 # coupling between Q and Bz

    np_gap_B = 100

    # coefficient to take into account additional multipliers in an oracle due to superpositions:
    # in 1-D fluid plasma X-wave problem: there is a superposition of states -> max multiplier is 2
    coef_superposition = 2 # might change for a different discretization scheme

    # --- prepare profiles and normalization ---
    dd_init = {
        "path-classical": path_root_classical,
        "path-qc": path_qc,
        "fname-start": fname_start,
        "r0": r0, "t0": t0,
        "n0": n0, "B0": B0, "R0": R0, "kx": kx,
        "nqx": nqx, "nx": nx, "nt": nt,
        "x_n0": x_n0, "x_nw": x_nw,
        "n0_2": n0_2, "x_n0_2": x_n0_2, "x_nw_2": x_nw_2,
        "B0_2": B0_2, "R0_2": R0_2, "x_B_left": x_B_left,
        "sel-norm": "c",
        "coef-coupling": coef_coupling,
        "coef-superposition": coef_superposition,
        "x-source": x0,
        "np-gap-B": np_gap_B
    }
    dd = fw.prepare_norm_profiles_outgoing(dd_init)
    return dd




def launch(dd):
    vars_data = fw.set_init_variables_Q_outgoing(dd["Q0"], dd["x-source"], dd)
    vars_data_res = fw.compute_Q_outgoing_DD(dd["coef-coupling"], vars_data, dd)
    vars_data_res["ids-Q"] = vars_data["ids-Q"]
    return vars_data_res







