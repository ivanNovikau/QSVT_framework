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


def choose_a_state(dd, choice, def_value=0):
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

    ch_state = choose_a_state(dd, choice)
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
    ch_state = choose_a_state(dd, vars_enc)

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


## Set a part of the "state" to a bit-array represented by an integer "int_repr".
# "state" is a list of size nq, where nq is a total number of qubits in a QSP circuit.
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
    ch_state = choose_a_state(dd, choice, -1)

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
# Remark: the first bit in a state has the highest priority.
def get_int_from_reg_state(dd, reg_name, input_state):
    n_reg = dd["regs"][reg_name]
    shift_reg = dd["reg-shifts"][reg_name]
    int_repr = 0
    for iq in range(n_reg):
        int_repr += 2**(n_reg-iq-1) * input_state[shift_reg + iq]
    return int_repr


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

        ampls_complex = self.get_complex(ampls)
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


    def get_complex(self, ampls):
        N = len(ampls)
        ampls_complex = np.zeros(N, dtype=np.complex)
        for i_state in range(N):
            one_ampl = ampls[i_state]
            ampls_complex[i_state] = np.complex(one_ampl["real"], one_ampl["imag"])
        return ampls_complex


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
        ch_state = choose_a_state(self.dd_, choice, -1)

        one_step_states = self.states_[id_t]["state"]
        one_step_ampls = self.states_[id_t]["ampls"]

        # normalize the state amplitudes:
        ampls_complex = self.get_complex(one_step_ampls)
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
        # prepare a preliminary state defined by "vars_enc":
        ch_state = self.choose_a_state(vars_enc)

        # create all possible combinations of the register "reg_x" in "ch_state":
        nx = self.dd_["regs"][reg_x]
        Nx = 2**nx
        ch_states_x = [None] * Nx
        for id_x in range(Nx):
            ch_states_x[id_x] = self.set_reg(ch_state, id_x, reg_x)
            # print(get_str_state(ch_states_x[id_x], dd["reg-nq"]))

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
            int_x = self.get_int_from_reg_state(reg_x, states_to_search[i_state])
            one_ampl = ampls_to_search[i_state]
            ampls[int_x] = np.complex(one_ampl["real"], one_ampl["imag"])
        return ampls


    def choose_a_state(self, choice, def_value=0):
        nq  = self.dd_["nq"]
        ch_state = [def_value] * nq
        for reg_name, reg_int in choice.items():
            bit_array = mix.find_bit_array_of_int(reg_int, self.dd_["regs"][reg_name])
            for i_bit in range(len(bit_array)):
                i_bit_pos = self.dd_["reg-shifts"][reg_name] + i_bit
                ch_state[i_bit_pos] = bit_array[i_bit]
        return ch_state


    ## Set a part of the "state" to a bit-array represented by the integer "int_repr".
    # "state" is a list of size nq, where nq is a total number of qubits in the QSP circuit.
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
        ch_state = choose_a_state(self.dd_, choice, -1)

        one_step_states = self.states_[id_t]["state"]
        one_step_ampls = self.states_[id_t]["ampls"]

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


    ## Return an integer representation of a state encoded in a register "reg_name"
    # in a given state "input_state".
    # Remark: the first bit in a state has the highest priority.
    def get_int_from_reg_state(self, reg_name, input_state):
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
            self.dd_["qsp"]["odd-angles"] = bg["odd-angles"][()]
            self.dd_["qsp"]["even-angles"] = bg["even-angles"][()]

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