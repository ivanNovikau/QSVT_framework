from asyncio import constants
from pickletools import read_unicodestringnl
import sys
import datetime
import re
from typing import no_type_check
import numpy as np
import time
import h5py

import os

import pylib.mix as mix
import pylib.Global_variables as GLO

def reload():
    mix.reload_module(mix)
    mix.reload_module(GLO)


def get_str_complex(v):
    ll = ""
    vr,  vi  = v["real"],  v["imag"]
    avr, avi = np.abs(vr), np.abs(vi)
    coef_zero = 1e-16

    if   vi > 0 and avi > coef_zero and avr > coef_zero:
        ll = "{:0.3e}+{:0.3e}j".format(vr, vi)
    elif vi < 0 and avi > coef_zero and avr > coef_zero:
        ll = "{:0.3e}{:0.3e}j".format(vr, vi)
    elif avi < coef_zero and avr > coef_zero:
        ll = "{:0.3e}".format(vr)
    elif avr < coef_zero and avi > coef_zero:
        ll = "{:0.3e}j".format(vi)
    else:
        ll = "0.0"
    return ll


def read_restart(dd):
    fname = dd["path"] + "/" + dd["fname"]
    print(f"Reading the file {fname}...")

    res_dd = {}
    with h5py.File(fname, "r") as f:
        fstate = f["state"]

        res_dd["ampls"] = np.array(fstate["real"]) + 1j*np.array(fstate["imag"])

    return res_dd


def open(dd):
    # dd is dictionary, must contain 
    # -> a path to hdf5 files dd["path"];
    # -> a name of a .hdf5 file to read dd["fname"];
    fname = dd["path"] + "/" + dd["fname"]

    print(f"Reading the file {fname}...")
    with h5py.File(fname, "r") as f:

        # --- basic data ---
        bg = f["basic"]

        dd["flag-restart"] = False

        dd["date-of-sim"]  = bg["date-of-simulation"][()].decode("utf-8")
        dd["project-name"] = bg["project-name"][()].decode("utf-8")
        dd["launch-path"] = bg["launch-path"][()].decode("utf-8")
        dd["path-inputs"] = bg["path-inputs"][()].decode("utf-8")

        dd["nq"] = bg["nq"][()]
        dd["na"] = bg["na"][()]

        reg_names = bg["register-names"][()].decode("utf-8").split(", ")
        reg_nq = bg["register-nq"][...]

        dd["reg-names"] = reg_names
        dd["reg-nq"] = reg_nq

        dd["regs"] = {}
        dd["reg-shifts"] = {}
        reg_shift = 0
        for i in range(len(reg_nq)):
            dd["regs"][reg_names[i]] = reg_nq[i]
            dd["reg-shifts"][reg_names[i]] = reg_shift
            reg_shift += reg_nq[i]

        t_grid = list(bg["time-grid"][...])
        dd["t"] = np.array([0] + t_grid) 

        dd["qsp"] = {}
        dd["qsp"]["dt"] = t_grid[1] - t_grid[0]
        dd["qsp"]["norm"] = bg["normalization-coef"][()]
        dd["qsp"]["prec"] = bg["qsp-initial-precision"][()]
        dd["qsp"]["angles"] = bg["qsp-angles-one-time-step"][()]

        # --- initial states ---
        st = f["states"]

        nt = np.size(dd["t"])
        dd["states"] = [None] * nt

        dd["states"][0] = {}
        dd["states"][0]["state"] = np.transpose(np.array(st["initial-states"])) 
        dd["states"][0]["ampls"] = np.array(st["initial-amplitudes"])

    print("Name of the simulation is", dd["project-name"])
    print("Simulation has been performed ", dd["date-of-sim"])

    return dd

## Open all output files from different restarts
def open_all(dd):
    # dd is dictionary, must contain 
    # -> a path to hdf5 files:        dd["path"];
    # -> a name of a project to read: dd["pname"];

    dd["flag-restart"] = True
    dd["t"] = np.zeros(1)

    id_res = -1
    while 1:
        id_res += 1

        strf = "start-" + str(id_res)
        dd[strf] = {}

        if id_res == 0:
            fname = dd["path"] + "/" + dd["pname"] + "_OUTPUT.hdf5"
        else:
            fname = dd["path"] + "/" + dd["pname"] + "_" + str(id_res) + "_OUTPUT.hdf5"

        if not os.path.isfile(fname):
            break

        print(f"\n--- Reading the file {fname}...")
        with h5py.File(fname, "r") as f:

            # --- basic data ---
            bg = f["basic"]

            dd[strf]["date-of-sim"]  = bg["date-of-simulation"][()].decode("utf-8")
            dd["project-name"] = bg["project-name"][()].decode("utf-8")
            dd[strf]["launch-path"] = bg["launch-path"][()].decode("utf-8")
            dd[strf]["path-inputs"] = bg["path-inputs"][()].decode("utf-8")

            dd["nq"] = bg["nq"][()]
            dd["na"] = bg["na"][()]

            reg_names = bg["register-names"][()].decode("utf-8").split(", ")
            reg_nq = bg["register-nq"][...]

            dd["reg-names"] = reg_names
            dd["reg-nq"] = reg_nq

            dd["regs"] = {}
            dd["reg-shifts"] = {}
            reg_shift = 0
            for i in range(len(reg_nq)):
                dd["regs"][reg_names[i]] = reg_nq[i]
                dd["reg-shifts"][reg_names[i]] = reg_shift
                reg_shift += reg_nq[i]

            # --- time grid ---
            t_grid = np.array(bg["time-grid"][...])
            t_grid += dd["t"][-1]
            dd["t"] = np.concatenate((dd["t"], t_grid))

            if id_res == 0:
                # --- QSP parameters ---
                dd["qsp"] = {}
                dd["qsp"]["dt"] = t_grid[1] - t_grid[0]
                dd["qsp"]["norm"] = bg["normalization-coef"][()]
                dd["qsp"]["prec"] = bg["qsp-initial-precision"][()]
                dd["qsp"]["angles"] = bg["qsp-angles-one-time-step"][()]

                # --- initial states ---
                st = f["states"]

                nt = np.size(dd["t"])
                dd["states"] = [None] * nt

                dd["states"][0] = {}
                dd["states"][0]["state"] = np.transpose(np.array(st["initial-states"])) 
                dd["states"][0]["ampls"] = np.array(st["initial-amplitudes"])
            else:
                dd["states"] += [None] * np.size(t_grid)

        print("Simulation has been performed ", dd[strf]["date-of-sim"])

    return dd


## Read all available quantum states and their amplitudes at every time step:
# results: 
# dd["states"][id_time_step]["state"] - 2-D np.array of states at a time step id_time_step;
# dd["states"][id_time_step]["ampls"] - 1-D np.array of amplitudes of the above states.
def read_all_output_states(dd):
    if not dd["flag-restart"]:
        fname = dd["path"] + "/" + dd["fname"]
        nt = np.size(dd["t"])
        with h5py.File(fname, "r") as f:
            for i in range(1,nt):
                str_step = "t-step-" + str(i-1)
                one_step_ampls = np.array(f["states/" + str_step + "--output-amplitudes"])
                one_step_states = np.array(f["states/" + str_step + "--output-states"])
                one_step_states = np.transpose(one_step_states)

                dd["states"][i] = {}
                dd["states"][i]["state"] = one_step_states
                dd["states"][i]["ampls"]  = one_step_ampls
    else:
        id_res = -1
        id_step_global = 0
        while 1:
            id_res += 1
            if id_res == 0:
                fname = dd["path"] + "/" + dd["pname"] + "_OUTPUT.hdf5"
            else:
                fname = dd["path"] + "/" + dd["pname"] + "_" + str(id_res) + "_OUTPUT.hdf5"
            if not os.path.isfile(fname):
                break

            print(f"\n--- Reading the file {fname}...")
            with h5py.File(fname, "r") as f:
                nt = np.size(f["basic/time-grid"][...])
                for i in range(nt):
                    id_step_global += 1
                    str_step = "t-step-" + str(i)
                    one_step_ampls = np.array(f["states/" + str_step + "--output-amplitudes"])
                    one_step_states = np.array(f["states/" + str_step + "--output-states"])
                    one_step_states = np.transpose(one_step_states)

                    dd["states"][id_step_global] = {}
                    dd["states"][id_step_global]["state"] = one_step_states
                    dd["states"][id_step_global]["ampls"] = one_step_ampls
    return


## Print a state and its amplitude:
#  E.g. dd_one_step can be dd["states"][0]
def print_state(dd, dd_one_step):
    state = dd_one_step["state"]
    ampls = dd_one_step["ampls"]

    nr, _ = state.shape

    print("Registers: ")
    print(dd["regs"])
    print("")
    for ir in range(nr):
        str_ampl = get_str_complex(ampls[ir])
        str_state = get_str_state(state[ir], dd["reg-nq"])
        print("{:>22s}   {:s}".format(str_ampl, str_state))
    return


def print_initial_states(dd):
    print("--- Initial state ---")
    print_state(dd, dd["states"][0])


def get_str_state(q, format_q):
    nq = np.size(q)

    ll = ""
    count_q = 0
    count_fq = 0
    for i in range(nq):
        if count_q == 0:
            ll += "|"
            
        count_q += 1
        ll += "{:d}".format(q[i])
        if count_q == format_q[count_fq]:
            ll += ">"
            count_fq += 1
            count_q = 0
    return ll


def create_mask(dd, choice, def_value=0):
    nq  = dd["nq"]
    ch_state = [def_value] * nq
    for reg_name, reg_int in choice.items():
        bit_array = mix.find_bit_array_of_int(reg_int, dd["regs"][reg_name])
        for i_bit in range(len(bit_array)):
            i_bit_pos = dd["reg-shifts"][reg_name] + i_bit
            ch_state[i_bit_pos] = bit_array[i_bit]
    return ch_state


## Return an amplitude and a state defined by the dictionary "choice":
# choice = {"reg_name_1": int_to_choose, "reg_name_2": int_to_choose, ...}.
# If a register is not defined in "choice", it is set to a zero state. 
# Return the chosen state and its amplitude at every time moment.
def get_state_on_t(dd, choice):
    print("Registers: ")
    print(dd["regs"])
    print("")

    ch_state = create_mask(dd, choice)
    print("Find amplitudes of the following state: ")
    print(get_str_state(ch_state, dd["reg-nq"]))

    nt = np.size(dd["t"])
    res_ampls = [np.complex(0.0)] * nt
    for i_step in range(nt):
         one_step_states = dd["states"][i_step]["state"]
         ns, _ = one_step_states.shape
         for i_state in range(ns):
             if np.array_equal(one_step_states[i_state], ch_state):
                 one_ampl = dd["states"][i_step]["ampls"][i_state]
                 res_ampls[i_step] = np.complex(one_ampl["real"], one_ampl["imag"])
    return res_ampls, ch_state


## Return a 1-D np.array (0.0, step, ... 1.0) defined by 
# the size (number of qubits) of the register with a name "reg_x".
def get_x_grid(dd, reg_x):
    reg_nq = dd["regs"][reg_x]
    N = 2**reg_nq
    x_grid = np.array(range(N))/(1.*(N-1))
    return x_grid


## Return amplitudes as a 2-D np.array (x,t) of a variable defined by "vars_enc".
# The space dependence on x is encoded in a register "reg_x".
# "vars_enc" is {"reg_name_1": int_to_choose, ...};
# "reg_x" is a name of a register where different combinations of qubits
#           correspond to different points on x.
# All other registers are set zero.
def get_amplitudes_xt(dd, vars_enc, reg_x):
    # prepare a preliminary state defined by "vars_enc":
    ch_state = create_mask(dd, vars_enc)

    # create all possible combinations of a register "reg_x" in "ch_state":
    Nt = np.size(dd["t"])
    nx = dd["regs"][reg_x]
    Nx = 2**nx
    ch_states_x = [None] * Nx
    for id_x in range(Nx):
        ch_states_x[id_x] = set_reg(dd, ch_state, id_x, reg_x)
        # print(get_str_state(ch_states_x[id_x], dd["reg-nq"]))

    # prepare a dictionary that defines a set of states to be considered at every time step:
    var_to_cons = {}
    for reg_name in dd["reg-names"]:
        if reg_name in vars_enc.keys():
            var_to_cons[reg_name] = vars_enc[reg_name]
            continue
        if reg_name != reg_x:
            var_to_cons[reg_name] = 0
            continue

    ampls = np.zeros((Nx, Nt), dtype=np.complex)
    for it in range(Nt):

        # consider only states for a chosen variable encoded by vars_enc:
        ampls_to_search, states_to_search = get_several_chosen_states_at_t1(dd, var_to_cons, it)

        # print(f"\nstates to consider at id_t = {it}:")
        # for one_state in states_to_search:
        #     print(get_str_state(one_state, dd["reg-nq"]))
        #     int_x = get_int_from_reg_state(dd, reg_x, one_state)
        #     print(f"int repres. is {int_x}\n")

        # every state in the considered set of states must correspond to one space point:
        nstates = len(states_to_search)
        for i_state in range(nstates):
            int_x = get_int_from_reg_state(dd, reg_x, states_to_search[i_state])
            one_ampl = ampls_to_search[i_state]
            ampls[int_x, it] = np.complex(one_ampl["real"], one_ampl["imag"])

            # print(f"at [ix,it] = [{int_x}, {it}]")
            # print(np.complex(one_ampl["real"], one_ampl["imag"]))
    return ampls


def get_total_prob(dd, id_t):
    ampls = dd["states"][id_t]["ampls"]

    tot_prob = 0.
    na = len(ampls)
    for i_state in range(na):
        one_ampl = ampls[i_state]
        one_ampl_compl = np.complex(one_ampl["real"], one_ampl["imag"])
        tot_prob += np.abs(one_ampl_compl)**2
    return tot_prob


## Set a part of the "state" to a bit-array represented by an integer "int_repr".
# "state" is a list of size nq, where nq is a total number of qubits in the circuit.
# "reg_name" defines a part of the array "state" that must be set to "int_repr".
# Return a new state as an 1-D np.array.
def set_reg(dd, state, int_repr, reg_name):
    n_reg = dd["regs"][reg_name]
    shift_reg = dd["reg-shifts"][reg_name]
    bit_array = mix.find_bit_array_of_int(int_repr, n_reg)

    res_state = np.array(state)
    for iq in range(n_reg):
        res_state[shift_reg + iq] = bit_array[iq]
    return res_state


## Return states and their amplitudes defined by the dictionary "choice":
# choice = {"reg_name_1": int_to_choose, "reg_name_2": int_to_choose, ...}.
# The dictionary indicates in which state must be a register of a circuit.
# If a register is not defined in the dictionary, then 
# the function returns all available states from this register.
# Return the chosen states and their amplitudes at a chosen time moment "id_t".
def get_several_chosen_states_at_t1(dd, choice, id_t):
    ch_state = create_mask(dd, choice, -1)

    one_step_states = dd["states"][id_t]["state"]
    one_step_ampls = dd["states"][id_t]["ampls"]

    nstates, _ = one_step_states.shape
    res_ampls = []
    res_states = []
    for i_state in range(nstates):
        flag_choose = True
        one_state = one_step_states[i_state]
        for i_qubit in range(dd["nq"]):
            if ch_state[i_qubit] == -1:
                continue
            if ch_state[i_qubit] != one_state[i_qubit]:
                flag_choose = False
                break
        if flag_choose:
            res_ampls.append(one_step_ampls[i_state])
            res_states.append(one_state)
    return res_ampls, res_states


## Return an integer representation of a state encoded in a register "reg_name"
# in a given state "input_state".
# Remark: the first bit in a state is the most significant.
# !!! Does not take into account the position of the register within the circuit !!!
# !!! Means, the same two bit-arrays, which lie in two different places of the circuit,
#            will return the same integer !!!
# !!! Yet, in general, more significant register must return a larger integers, just because
#   the register has more significant qubits !!!
def get_int_from_reg_state(dd, reg_name, input_state):
    n_reg = dd["regs"][reg_name]
    shift_reg = dd["reg-shifts"][reg_name]
    int_repr = 0
    for iq in range(n_reg):
        int_repr += 2**(n_reg-iq-1) * input_state[shift_reg + iq]
    return int_repr


def get_complex(ampls):
    N = len(ampls)
    ampls_complex = np.zeros(N, dtype=np.complex)
    for i_state in range(N):
        one_ampl = ampls[i_state]
        ampls_complex[i_state] = np.complex(one_ampl["real"], one_ampl["imag"])
    return ampls_complex


def compare_states(output_state_dict, qiskit_ampls, qiskit_states, err_comp = 1e-6, flag_print = False):
    states = output_state_dict["state"]
    ampls = get_complex(output_state_dict["ampls"])

    n_states_, _ = states.shape
    if not (n_states_ == len(qiskit_states)):
        print("different number of states")
        return

    # compare state by state:
    for i_state in range(n_states_):
        ampl_one = ampls[i_state]
        diff_real = ampl_one.real - qiskit_ampls[i_state].real
        diff_imag = ampl_one.imag - qiskit_ampls[i_state].imag

        if flag_print:
            print("\n-----------------------------")
            print(list(states[i_state]))
            print(qiskit_states[i_state])
        if not (states[i_state] == qiskit_states[i_state]).all():
            print("the states are different")
            return

        if flag_print:
            print(diff_real)
            print(diff_imag)
        if np.abs(diff_real) < err_comp and np.abs(diff_imag) < err_comp:
            continue
        else:
            print("the amplitudes are different")
            return
    print("the states are close to each other with the imposed error.")
    return


# calculate the total probability of the state defined only by the list_qubits;
def calc_tot_prob_wrt(state_dict, list_qubits):
    ampls = get_complex(state_dict["ampls"])
    n_states, _ = state_dict["state"].shape
    new_state_dict = {"states": [], "probs": []}
    states_checked = []
    for i_state in range(n_states):
        prob_state = 0
        one_state = state_dict["state"][i_state][list_qubits]
        if list(one_state) in states_checked:
            continue
        for j_state in range(n_states):
            if(one_state == state_dict["state"][j_state][list_qubits]).all():
                prob_state += np.abs(ampls[j_state])**2
        new_state_dict["states"].append(one_state)
        new_state_dict["probs"].append(prob_state)
        states_checked.append(list(one_state))
    new_state_dict["states"] = np.array(new_state_dict["states"])
    new_state_dict["probs"]  = np.array(new_state_dict["probs"])
    return new_state_dict


# probs_int - data structure obtained from calc_tot_prob_wrt:
def pe_form_hist(probs_int):
    n_states, n_qubits = probs_int["states"].shape
    hist_probs = np.zeros(n_states)
    hist_phases = np.zeros(n_states)
    count_hist = -1
    for i_state in range(n_states):
        count_hist += 1
        cl_state = probs_int["states"][i_state]
        int_1 = mix.find_int_from_bit_array(cl_state)
        hist_phases[count_hist] = 2*np.pi * int_1 / 2**n_qubits
        hist_probs[count_hist]  = probs_int["probs"][i_state]
    hist_res = {
        "probs": np.array(hist_probs),
        "phases": np.array(hist_phases)
    }
    return hist_res


# -------------------------------------------------------------------------------
# --- Basic class for measurements ---
# -------------------------------------------------------------------------------
class Meas__:
    # project name:
    pname_ = ""

    # path to the project:
    path_ = ""

    # information about the project:
    dd_ = None

    # states:  
    states_ = None


    def __init__(self):
        self.dd_ = {}

        # initial and final states:
        
        return 


    def open(self):
        fname = self.path_ + "/" + self.pname_ + "_OUTPUT.hdf5"
        self.dd_["fname"] = fname

        print(f"Reading the file {fname}...")
        with h5py.File(fname, "r") as f:

            # --- basic data ---
            bg = f["basic"]

            self.dd_["flag-restart"] = False

            self.dd_["date-of-sim"]  = bg["date-of-simulation"][()].decode("utf-8")
            self.dd_["project-name"] = bg["project-name"][()].decode("utf-8")
            self.dd_["launch-path"] = bg["launch-path"][()].decode("utf-8")
            self.dd_["path-inputs"] = bg["path-inputs"][()].decode("utf-8")

            self.dd_["nq"] = bg["nq"][()]
            self.dd_["na"] = bg["na"][()]

            reg_names = bg["register-names"][()].decode("utf-8").split(", ")
            reg_nq = bg["register-nq"][...]

            self.dd_["reg-names"] = reg_names
            self.dd_["reg-nq"] = reg_nq

            self.dd_["regs"] = {}
            self.dd_["reg-shifts"] = {}
            reg_shift = 0
            for i in range(len(reg_nq)):
                self.dd_["regs"][reg_names[i]] = reg_nq[i]
                self.dd_["reg-shifts"][reg_names[i]] = reg_shift
                reg_shift += reg_nq[i]

            self.dd_["qsp"] = {}
            self.dd_["qsp"]["norm"] = bg["normalization-coef"][()]
            self.dd_["qsp"]["prec"] = bg["qsvt-error"][()]

            # --- initial states ---
            st = f["states"]
            self.states_[0]["state"] = np.transpose(np.array(st["initial-states"])) 
            self.states_[0]["ampls"] = np.array(st["initial-amplitudes"])

        print("Name of the simulation is", self.dd_["project-name"])
        print("Simulation has been performed ", self.dd_["date-of-sim"])
        return


    def print_state(self, id_step):
        # id_step = 0 - initial state
        # id_step = -1 - last state
        state = self.states_[id_step]["state"]
        ampls = self.states_[id_step]["ampls"]

        nr, _ = state.shape

        print("Registers: ")
        print(self.dd_["regs"])
        print("")
        for ir in range(nr):
            str_ampl = get_str_complex(ampls[ir])
            str_state = get_str_state(state[ir], self.dd_["reg-nq"])
            print("{:>22s}   {:s}".format(str_ampl, str_state))
        return


    def calc_rel_state(self, id_step):
        # id_step = 0 - initial state
        # id_step = -1 - last state
        state = self.states_[id_step]["state"]
        ampls = self.states_[id_step]["ampls"]

        ampls_complex = get_complex(ampls)
        norm_ampl = np.sqrt(np.sum(np.abs(ampls_complex)**2))
        ampls_norm = ampls_complex / norm_ampl
        return ampls_norm, state


    def print_rel_state(self, id_step, flag_abs=False):
        ampls_norm, state = self.calc_rel_state(id_step)
        nr, _ = state.shape

        print("Registers: ")
        print(self.dd_["regs"])
        print("")
        for ir in range(nr):
            str_state = get_str_state(state[ir], self.dd_["reg-nq"])
            if flag_abs:
                print("{:40.3e}   {:s}".format(np.abs(ampls_norm[ir]), str_state))
            else:
                print("{:40.3e}   {:s}".format(ampls_norm[ir], str_state))
        return


    ## 
    # -> choice = {"reg_name_1": int_1, "reg_name_2": int_2, ...}.
    # One takes all states, which satisfy:
    # reg_1 = |int_1>, reg_2 = |int_2>, ...
    # If the register is not defined in the dictionary, then 
    # the function returns all available states from this register.
    # -> Return the chosen states and their normalized amplitudes at the chosen time moment id_t.
    # Normalized (relative) means each amplitude is normalized to the total amplitude of all
    # available states. 
    # Usually, the available states are the states with zero ancillae registers.
    def get_rel_ampl(self, choice, id_t):
        ch_state = self.create_mask(choice, -1)

        one_step_states = self.states_[id_t]["state"]
        one_step_ampls = self.states_[id_t]["ampls"]

        # normalize the state amplitudes:
        ampls_complex = get_complex(one_step_ampls)
        norm_ampl = np.sqrt(np.sum(np.abs(ampls_complex)**2))
        one_step_ampls = ampls_complex / norm_ampl

        # choose the necessary states:
        nstates, _ = one_step_states.shape
        res_ampls = []
        res_states = []
        for i_state in range(nstates):
            flag_choose = True
            one_state = one_step_states[i_state]
            for i_qubit in range(self.dd_["nq"]):
                if ch_state[i_qubit] == -1:
                    continue
                if ch_state[i_qubit] != one_state[i_qubit]:
                    flag_choose = False
                    break
            if flag_choose:
                res_ampls.append(one_step_ampls[i_state])
                res_states.append(one_state)
        return res_ampls, res_states


    ## Return a 1-D np.array (0.0, step, ... 1.0) defined by 
    # the number of qubits in the register with the name "reg_x".
    def get_x_grid(self, reg_x):
        reg_nq = self.dd_["regs"][reg_x]
        N = 2**reg_nq
        x_grid = np.array(range(N))/(1.*(N-1))
        return x_grid


    ## Return the variable defined by "vars_enc" as a function of x.
    # The space dependence on x is encoded in the register "reg_x".
    # "vars_enc" is {"reg_name_1": int_to_choose, ...};
    # "reg_x" is the register name  where different combinations of qubits
    #           correspond to different points on x.
    # All other registers are set to zero.
    def get_var_x(self, vars_enc, reg_x):
        # # prepare a preliminary state defined by "vars_enc":
        # ch_state = self.create_mask(vars_enc)

        # # create all possible combinations of the register "reg_x" in "ch_state":
        nx = self.dd_["regs"][reg_x]
        Nx = 2**nx
        # ch_states_x = [None] * Nx
        # for id_x in range(Nx):
        #     ch_states_x[id_x] = self.set_reg(ch_state, id_x, reg_x)
        #     # print(get_str_state(ch_states_x[id_x], dd["reg-nq"]))

        # prepare a dictionary that defines a set of states to be considered at every time step:
        var_to_cons = {}
        for reg_name in self.dd_["reg-names"]:
            if reg_name in vars_enc.keys():
                var_to_cons[reg_name] = vars_enc[reg_name]
                continue
            if reg_name != reg_x:
                var_to_cons[reg_name] = 0
                continue

        ampls = np.zeros(Nx, dtype=np.complex)

        # consider only states for the chosen variable encoded by vars_enc:
        ampls_to_search, states_to_search = self.get_several_chosen_states_at_t1(var_to_cons, -1)

        # every state in the considered set of states must correspond to one space point:
        nstates = len(states_to_search)
        for i_state in range(nstates):
            int_x = self.convert_reg_state_to_int(reg_x, states_to_search[i_state])
            one_ampl = ampls_to_search[i_state]
            ampls[int_x] = np.complex(one_ampl["real"], one_ampl["imag"])
        return ampls

    ## Return the variable defined by "vars_enc" as a function of x1 and x2.
    # The dependence on x1 (x2) is encoded in the register reg_x1 (reg_x2).
    # "vars_enc" is {"reg_name_1": int_to_choose, ...};
    # "reg_x1" (reg_x2) is the register name where different combinations of qubits
    #           correspond to different points on x1 (x2).
    # All other registers are set to zero.
    def get_var_x1x2(self, vars_enc, reg_x1, reg_x2):
        Nx1 = 2**self.dd_["regs"][reg_x1]
        Nx2 = 2**self.dd_["regs"][reg_x2]

        # prepare a dictionary that defines a set of states to be considered at every time step:
        var_to_cons = {}
        for reg_name in self.dd_["reg-names"]:
            if reg_name in vars_enc.keys():
                var_to_cons[reg_name] = vars_enc[reg_name]
                continue
            if reg_name != reg_x1 and reg_name != reg_x2:
                var_to_cons[reg_name] = 0
                continue
        ampls = np.zeros((Nx1, Nx2), dtype=np.complex)

        # consider only states for the chosen variable encoded by vars_enc:
        ampls_to_search, states_to_search = self.get_several_chosen_states_at_t1(var_to_cons, -1)

        # every state in the considered set of states must correspond to one space point:
        nstates = len(states_to_search)
        for i_state in range(nstates):
            int_x1 = self.convert_reg_state_to_int(reg_x1, states_to_search[i_state])
            int_x2 = self.convert_reg_state_to_int(reg_x2, states_to_search[i_state])
            one_ampl = ampls_to_search[i_state]
            ampls[int_x1, int_x2] = np.complex(one_ampl["real"], one_ampl["imag"])
        return ampls


    ## Create an 1D-array of size 2**nq (nq - number of qubits in the circuit),
    # where corresponding pieces of the array are filled by binary representations of
    # the register states from the "choice" map.
    # The rest of the array elements are filled with the "def_value". 
    def create_mask(self, choice, def_value=0):
        nq  = self.dd_["nq"]
        ch_state = [def_value] * nq
        for reg_name, reg_int in choice.items():
            bit_array = mix.find_bit_array_of_int(reg_int, self.dd_["regs"][reg_name])
            for i_bit in range(len(bit_array)):
                i_bit_pos = self.dd_["reg-shifts"][reg_name] + i_bit
                ch_state[i_bit_pos] = bit_array[i_bit]
        return ch_state


    ## Set a part of the "state" to a bit-array represented by the integer "int_repr".
    # "state" is a list of size nq, where nq is a total number of qubits in the circuit.
    # "reg_name" defines a part of the array "state" that must be set to "int_repr".
    # Return a new state as an 1-D np.array.
    def set_reg(self, state, int_repr, reg_name):
        n_reg = self.dd_["regs"][reg_name]
        shift_reg = self.dd_["reg-shifts"][reg_name]
        bit_array = mix.find_bit_array_of_int(int_repr, n_reg)

        res_state = np.array(state)
        for iq in range(n_reg):
            res_state[shift_reg + iq] = bit_array[iq]
        return res_state


    ## Return states and their amplitudes defined by the dictionary "choice":
    # choice = {"reg_name_1": int_to_choose, "reg_name_2": int_to_choose, ...}.
    # The dictionary indicates which state the register must have.
    # If the register is not defined in the dictionary, then 
    # the function returns all available states from this register.
    # Return the chosen states and their amplitudes at the chosen time moment "id_t".
    def get_several_chosen_states_at_t1(self, choice, id_t):
        ch_state = self.create_mask(choice, -1)

        one_step_states = self.states_[id_t]["state"]
        one_step_ampls  = self.states_[id_t]["ampls"]

        nstates, _ = one_step_states.shape
        res_ampls = []
        res_states = []
        for i_state in range(nstates):
            flag_choose = True
            one_state = one_step_states[i_state]
            for i_qubit in range(self.dd_["nq"]):
                if ch_state[i_qubit] == -1:
                    continue
                if ch_state[i_qubit] != one_state[i_qubit]:
                    flag_choose = False
                    break
            if flag_choose:
                res_ampls.append(one_step_ampls[i_state])
                res_states.append(one_state)
        return res_ampls, res_states


    ## The "input_state" is 1-D array of size 2**nq (nq - number of qubits in the circuit);
    # The register name "reg_name" chooses a piece of the array "input_state".
    # This piece contains a bit-array.
    # Taking into account nq and (but the position of the bit-array within the array "input_state"),
    # the function returns the integer representation of the bit-array.
    # THe zeroth bit in the register is assumed to be the most significant.
    def convert_reg_state_to_int(self, reg_name, input_state):
        n_reg = self.dd_["regs"][reg_name]
        shift_reg = self.dd_["reg-shifts"][reg_name]
        int_repr = 0
        for iq in range(n_reg):
            int_repr += 2**(n_reg-iq-1) * input_state[shift_reg + iq]
        return int_repr


# -------------------------------------------------------------------------------
# --- Matrix inversion ---
# -------------------------------------------------------------------------------
class MeasInverse__(Meas__):
    def __init__(self):
        Meas__.__init__(self)
        self.states_ = [{}, {}]
        return 


    def open(self):
        Meas__.open(self)
        with h5py.File(self.dd_["fname"], "r") as f:
            bg = f["basic"]

            self.dd_["kappa"] = bg["kappa"][()]
            self.dd_["eps"] = bg["qsvt-error"][()]
            self.dd_["qsp"]["odd-angles"] = bg["odd-angles"][()]
        return


    def read_output_states(self):
        with h5py.File(self.dd_["fname"], "r") as f:
            st = f["states"]
            self.states_[-1]["state"] = np.transpose(np.array(st["output-states"])) 
            self.states_[-1]["ampls"] = np.array(st["output-amplitudes"])
        return


# -------------------------------------------------------------------------------
# --- Time Evolution ---
# -------------------------------------------------------------------------------
class MeasDyn__(Meas__):
    def __init__(self):
        Meas__.__init__(self)
        self.states_ = [{}]
        return 


    def open(self):
        Meas__.open(self)
        with h5py.File(self.dd_["fname"], "r") as f:
            bg = f["basic"]

            t_grid = list(bg["time-grid"][...])
            self.dd_["t"] = np.array([0] + t_grid) 

            self.dd_["qsp"]["dt"] = bg["dt"][()]
            # self.dd_["qsp"]["odd-angles"] = bg["odd-angles"][()]
            # self.dd_["qsp"]["even-angles"] = bg["even-angles"][()]

            self.states_ = self.states_ + [{}] * len(t_grid)
        return

    
    ## Read all available quantum states and their amplitudes at every time step:
    # results: 
    # self.dd_["states"][id_time_step]["state"] - 2-D np.array of states at the time step id_time_step;
    # self.dd_["states"][id_time_step]["ampls"] - 1-D np.array of amplitudes of the above states.
    def read_output_states(self):
        if not self.dd_["flag-restart"]:
            fname = self.dd_["fname"]
            nt = np.size(self.dd_["t"])
            with h5py.File(fname, "r") as f:
                for i in range(1,nt):
                    str_step = "t-step-" + str(i-1)
                    one_step_ampls = np.array(f["states/" + str_step + "--output-amplitudes"])
                    one_step_states = np.array(f["states/" + str_step + "--output-states"])
                    one_step_states = np.transpose(one_step_states)

                    self.states_[i] = {}
                    self.states_[i]["state"] = one_step_states
                    self.states_[i]["ampls"]  = one_step_ampls
        else:
            id_res = -1
            id_step_global = 0
            while 1:
                id_res += 1
                if id_res == 0:
                    fname = self.dd_["fname"]
                else:
                    fname = self.path_ + "/" + self.pname_ + "_" + str(id_res) + "_OUTPUT.hdf5"
                if not os.path.isfile(fname):
                    break

                print(f"\n--- Reading the file {fname}...")
                with h5py.File(fname, "r") as f:
                    nt = np.size(f["basic/time-grid"][...])
                    for i in range(nt):
                        id_step_global += 1
                        str_step = "t-step-" + str(i)
                        one_step_ampls = np.array(f["states/" + str_step + "--output-amplitudes"])
                        one_step_states = np.array(f["states/" + str_step + "--output-states"])
                        one_step_states = np.transpose(one_step_states)

                        self.states_[id_step_global] = {}
                        self.states_[id_step_global]["state"] = one_step_states
                        self.states_[id_step_global]["ampls"] = one_step_ampls
        return


# -------------------------------------------------------------------------------
# --- Read output data from oracle ---
# -------------------------------------------------------------------------------
class MeasOracle__:
    # project name:
    pname_ = ""

    # path to the project:
    path_ = ""

    # information about the project:
    dd_ = None

    # constants:
    constants_ = None

    # states:  
    n_init_states_ = None
    init_states_ = None
    output_all_states_ = None
    output_zero_anc_states_ = None


    def __init__(self):
        self.dd_ = {}
        return 


    def open(self):
        fname = self.path_ + "/" + self.pname_ + "_circuit_OUTPUT.hdf5"
        self.dd_["fname"] = fname

        print(f"Reading the file {fname}...")
        with h5py.File(fname, "r") as f:

            # --- basic data ---
            bg = f["basic"]

            self.dd_["date-of-sim"]  = bg["date-of-simulation"][()].decode("utf-8")
            self.dd_["project-name"] = bg["project-name"][()].decode("utf-8")
            self.dd_["launch-path"] = bg["launch-path"][()].decode("utf-8")
            self.dd_["path-inputs"] = bg["path-inputs"][()].decode("utf-8")

            self.dd_["nq"] = bg["nq"][()]
            self.dd_["na"] = bg["na"][()]

            reg_names = bg["register-names"][()].decode("utf-8").split(", ")
            reg_nq = bg["register-nq"][...]

            self.dd_["reg-names"] = reg_names
            self.dd_["reg-nq"] = reg_nq

            self.dd_["regs"] = {}
            self.dd_["reg-shifts"] = {}
            reg_shift = 0
            for i in range(len(reg_nq)):
                self.dd_["regs"][reg_names[i]] = reg_nq[i]
                self.dd_["reg-shifts"][reg_names[i]] = reg_shift
                reg_shift += reg_nq[i]

            # --- constants ---
            self.constants_ = {}
            for field in f["constants"]:
                self.constants_[field] = f["constants"][field][()]

            # --- initial states ---
            st = f["states"]
            self.n_init_states_ = st["n-init-states"][()]
            self.init_states_            = [{}] * self.n_init_states_
            self.output_all_states_      = [{}] * self.n_init_states_
            self.output_zero_anc_states_ = [{}] * self.n_init_states_
            for ii in range(self.n_init_states_):
                self.init_states_[ii]["state"] = np.transpose(np.array(st["initial-states-{:d}".format(ii)])) 
                self.init_states_[ii]["ampls"] = np.array(st["initial-amplitudes-{:d}".format(ii)])
                self.output_all_states_[ii]["state"] = np.transpose(np.array(st["output-all-states-{:d}".format(ii)])) 
                self.output_all_states_[ii]["ampls"] = np.array(st["output-all-amplitudes-{:d}".format(ii)])

                line_state = "output-zero-anc-states-{:d}".format(ii)
                if line_state in st:
                    self.output_zero_anc_states_[ii]["state"] = np.transpose(np.array(st[line_state])) 
                    self.output_zero_anc_states_[ii]["ampls"] = np.array(st["output-zero-anc-amplitudes-{:d}".format(ii)])

        print("Name of the simulation is", self.dd_["project-name"])
        print("Simulation has been performed ", self.dd_["date-of-sim"])
        return


    def print_full_states(self):
        print("Number of initial states: {:d}".format(self.n_init_states_))

        print("\nRegisters: ")
        print(self.dd_["regs"])
        print()

        for ii in range(self.n_init_states_):
            print("\n-------------------------------------")
            print("--- Initial state: {:d}".format(ii))
            state = self.init_states_[ii]["state"]
            ampls = self.init_states_[ii]["ampls"]
            nr, _ = state.shape
            for ir in range(nr):
                str_ampl = get_str_complex(ampls[ir])
                str_state = get_str_state(state[ir], self.dd_["reg-nq"])
                print("{:>22s}   {:s}".format(str_ampl, str_state))

            print("\n -- full output state --")
            state = self.output_all_states_[ii]["state"]
            ampls = self.output_all_states_[ii]["ampls"]
            nr, _ = state.shape
            for ir in range(nr):
                str_ampl = get_str_complex(ampls[ir])
                str_state = get_str_state(state[ir], self.dd_["reg-nq"])
                print("{:>22s}   {:s}".format(str_ampl, str_state))
        return


    def print_zero_anc_states(self):
        print("Number of initial states: {:d}".format(self.n_init_states_))

        print("\nRegisters: ")
        print(self.dd_["regs"])
        print()

        for ii in range(self.n_init_states_):
            print("\n-------------------------------------")
            print("--- Initial state: {:d}".format(ii))
            state = self.init_states_[ii]["state"]
            ampls = self.init_states_[ii]["ampls"]
            nr, _ = state.shape
            for ir in range(nr):
                str_ampl = get_str_complex(ampls[ir])
                str_state = get_str_state(state[ir], self.dd_["reg-nq"])
                print("{:>22s}   {:s}".format(str_ampl, str_state))

            print("\n -- zero-ancilla output state --")
            if(self.output_zero_anc_states_[ii]):
                state = self.output_zero_anc_states_[ii]["state"]
                ampls = self.output_zero_anc_states_[ii]["ampls"]
                nr, _ = state.shape
                for ir in range(nr):
                    str_ampl = get_str_complex(ampls[ir])
                    str_state = get_str_state(state[ir], self.dd_["reg-nq"])
                    print("{:>22s}   {:s}".format(str_ampl, str_state))
        return


    


    