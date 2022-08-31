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

    probs_ = None
    qubits_probs_ = None

    # states and amplitudes for the given input state;
    work_ampls_ = None
    work_states_ = None  


    def __init__(self):
        self.dd_ = {}
        return 
        

    def open(self):
        fname = self.path_ + "/" + self.pname_ + "_OUTPUT.hdf5"
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

            # --- initial and output states ---
            st = f["states"]
            self.n_init_states_ = st["n-init-states"][()]
            self.init_states_            = [{}] * self.n_init_states_
            self.output_all_states_      = [{}] * self.n_init_states_
            self.output_zero_anc_states_ = [{}] * self.n_init_states_
            for ii in range(self.n_init_states_):
                self.init_states_[ii]["state"] = np.transpose(np.array(st["initial-states-{:d}".format(ii)])) 
                self.init_states_[ii]["ampls"] = np.array(st["initial-amplitudes-{:d}".format(ii)])

                line_state = "output-all-states-{:d}".format(ii)
                if line_state in st:
                    self.output_all_states_[ii]["state"] = np.transpose(np.array(st["output-all-states-{:d}".format(ii)])) 
                    self.output_all_states_[ii]["ampls"] = np.array(st["output-all-amplitudes-{:d}".format(ii)])

                line_state = "output-zero-anc-states-{:d}".format(ii)
                if line_state in st:
                    self.output_zero_anc_states_[ii]["state"] = np.transpose(np.array(st[line_state])) 
                    self.output_zero_anc_states_[ii]["ampls"] = np.array(st["output-zero-anc-amplitudes-{:d}".format(ii)])

            # --- probabilities ---
            line_group = "probabilities"
            if line_group in f:
                pr = f[line_group]
                self.probs_        = np.array(pr["probs"])
                self.qubits_probs_ = np.array(pr["qubits"])

        print("Name of the simulation is", self.dd_["project-name"])
        print("Simulation has been performed ", self.dd_["date-of-sim"])
        return


    def set_zero_ancillae_work_states(self, id_input_state = 0):
        self.work_states_ = self.output_zero_anc_states_[id_input_state]["state"]
        self.work_ampls_  = self.output_zero_anc_states_[id_input_state]["ampls"]
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


    def read_qsvt(self):
        with h5py.File(self.dd_["fname"], "r") as f:
            gr = f["qsvt"]

            self.dd_["qsvt"] = {}
            self.dd_["qsvt"]["type"] = gr["type"][()].decode("utf-8")
            self.dd_["qsvt"]["eps"] = gr["eps"][()]
            if self.dd_["qsvt"]["type"] == "matrix-inversion":
                self.read_qsvt_matrix_inversion(gr)
            if self.dd_["qsvt"]["type"] == "hamiltonian-sim":
                self.read_qsvt_hamiltonian_sim(gr)   
        return


    def read_qsvt_matrix_inversion(self, gr):
        self.dd_["qsvt"]["kappa"] = gr["kappa"][()]
        print("kappa: {:0.3f}".format(self.dd_["qsvt"]["kappa"]))
        return


    def read_qsvt_hamiltonian_sim(self, gr):
        self.dd_["qsvt"]["dt"] = gr["dt"][()]
        self.dd_["qsvt"]["nt"] = gr["nt"][()]
        print("dt: {:0.3f}".format(self.dd_["qsvt"]["dt"]))
        print("nt: {:0.3f}".format(self.dd_["qsvt"]["nt"]))
        return


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
        nx = self.dd_["regs"][reg_x]
        Nx = 2**nx
        ampls = np.zeros(Nx, dtype=np.complex)

        # prepare a dictionary that defines a set of states to be considered:
        var_to_cons = {}
        for reg_name in self.dd_["reg-names"]:
            if reg_name in vars_enc.keys():
                var_to_cons[reg_name] = vars_enc[reg_name]
                continue
            if reg_name != reg_x:
                var_to_cons[reg_name] = 0
                continue

        # consider only states for the chosen variable encoded by vars_enc:
        ampls_to_search, states_to_search = self.get_several_chosen_work_states(var_to_cons)

        # every state in the considered set of states must correspond to one space point:
        nstates = len(states_to_search)
        for i_state in range(nstates):
            int_x = self.convert_reg_state_to_int(reg_x, states_to_search[i_state])
            one_ampl     = ampls_to_search[i_state]
            ampls[int_x] = np.complex(one_ampl["real"], one_ampl["imag"])
        return ampls


    ## Return states and their amplitudes defined by the dictionary "choice":
    # choice = {"reg_name_1": int_to_choose, "reg_name_2": int_to_choose, ...}.
    # The dictionary indicates which state the register must have.
    # If the register is not defined in the dictionary, then 
    # the function returns all available states from this register.
    def get_several_chosen_work_states(self, choice):
        ch_state = self.create_mask(choice, -1)

        one_step_states = self.work_states_
        one_step_ampls  = self.work_ampls_

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


    ## The "input_state" is 1-D array of size 2**nq (nq - number of qubits in the circuit);
    # The register name "reg_name" chooses a piece of the array "input_state".
    # This piece contains a bit-array.
    # The function returns the integer representation of the bit-array.
    # The zero-th bit in the register is assumed to be the most significant.
    def convert_reg_state_to_int(self, reg_name, input_state):
        n_reg = self.dd_["regs"][reg_name]
        shift_reg = self.dd_["reg-shifts"][reg_name]
        int_repr = 0
        for iq in range(n_reg):
            int_repr += 2**(n_reg-iq-1) * input_state[shift_reg + iq]
        return int_repr