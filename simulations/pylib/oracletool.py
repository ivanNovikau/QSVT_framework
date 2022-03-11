import tkinter as tk
import sys
import datetime
import re
import qiskit
from qiskit.tools.visualization import circuit_drawer
import numpy as np
import time

import subprocess

try:
    import pylib.mix as mix
    import pylib.Global_variables as GLO
    import pylib.draw_circuit as dcm
    # import pylib.OTFrames as otf
except:
    import mix as mix
    import Global_variables as GLO
    import draw_circuit as dcm
    # import OTFrames as otf

class OracleTool:
    timer_start_ = None

    name_ = "" # project name
    path_ = "./" # path to an input file
    input_file_ext_ = ".oracle"
    output_py_file_ext_ = ".pylog"
    circuit_ext_ = ".circuit"

    log_name_ = "" # name of an output (log) file

    keyword_comment_ = "//"

    fig_circuit_extention_ = ".png"

    

    # path to the c++ executable file oracletool
    path_c_oracletool_ = ""

    # number of MPI processors:
    n_mpi_procs_ = 1

    # --- Dictionary with Circuits ---
    # {
    #   name_of_circuit: {
    #       "name"        - circuit name
    #       "nq"          - number of qubits in the circuit
    #       "n_regs"      - number of registers
    #       "obj"         - QuantumCircuit 
    #       "flag_draw"   - if True, draw the circuit
    #       "flag_redraw" - if True, redraw after launching C++ QPlasma  
    #       "regs": {
    #           name_of_a_reg: {
    #               "name" - name of the register
    #               "obj" - QuantumRegister
    #               "nq" - number of qubits in the register
    #               "id_first_qubit" - global position of the first qubit of the register 
    #               "location" - position of the register from the top: 0 - at the top, 1 - next register, etc.
    #           }    
    #       }
    #   }
    # }
    circuits_ = None 

    # --- Constants ---
    constants_ = None

    # how to launch the program: draw, init_state, full
    sel_launch_ = "full"
    available_launch_selectors_ = ["draw", "init_state", "full", "full_iterator"]
    
    # -----------------------------------------------------------
    # --- CONSTRUCTOR ---
    # -----------------------------------------------------------
    # python3 
    #   [path_to_oracletool.py]/oracletool.py 
    #   [project_name] 
    #   [path_to_c_exe] 
    #   [working_directory_with_an_input_file]
    #   [n_mpi_processors]
    def __init__(self):
        self.timer_start_ = time.perf_counter()

        self.circuits_ = {}
        self.constants_ = {}

        n_args = len(sys.argv)

        # read project name
        if n_args < 2:
            print("Error: no project name is given as an argument to the program.")
            sys.exit(-1)
        self.name_ = sys.argv[1]
        print(f"Project name: {self.name_}")

        # read path to a c++ oracletool executable file:
        if n_args < 3:
            print("Error: no path to a c++ oracletool executable file is given as an argument to the program.")
            sys.exit(-1)
        self.path_c_oracletool_ = sys.argv[2]
        print(f"Path to an executable QPlasma file: {self.path_c_oracletool_}")

        # read path to an input file
        if n_args < 4:
            print("Error: no path to an input file is given as an argument to the program.")
            sys.exit(-1)
        self.path_ = sys.argv[3]
        print(f"Path to an input file: {self.path_}")

        # number of MPI processors:
        if n_args > 4:
            self.n_mpi_procs_ = int(sys.argv[4])
        print(f"Number of MPI processors: {self.n_mpi_procs_}")

        # create a log file:
        self.log_name_ = self.path_ + "/" + self.name_ + self.output_py_file_ext_
        ff = open(self.log_name_, "w")
        ff.write(datetime.datetime.now().strftime("%d-%b-%Y (%H:%M:%S.%f)") + "\n")
        ff.write("Project name: {}\n".format(self.name_))
        ff.write("Path to an input file: {}\n".format(self.path_))
        ff.close()
  
    # -----------------------------------------------------------
    # --- CLOSE ORACLETOOL ---
    # -----------------------------------------------------------
    def close_oracletool(self):
        self.print_log("\n*** Close the OracleTool ***")
        temp  = time.perf_counter()
        self.print_log(f"Time of simulation: {temp - self.timer_start_:0.4f} seconds")
        return
        
    # -----------------------------------------------------------
    # --- READ INPUT FILE ---
    # -----------------------------------------------------------
    def read_input_file(self):
        self.print_log("\n*** Reading an input file... ***")
        fname = self.path_ + "/" + self.name_ + self.input_file_ext_
        try:
            ff = open(fname, "r")
        except:
            raise Exception(f"File {fname} is not found.")

        data = ff.read()
        ff.close()
        data_words = self.remove_comments_empty_lines(data)

        id_word = 0
        n_words = len(data_words)
        while id_word < n_words:
            word = data_words[id_word]

            if mix.compare_two_strings(word, "CONSTANTS"):
                self.print_log("\nReading constants...")
                id_word += 1
                id_word = self.read_constants(data_words, id_word)
                if id_word is None:
                    raise Exception("End of the CONSTANTS section is not found.")
                continue

            if mix.compare_two_strings(word, "CIRCUITS_DECLARATION"):
                self.print_log("\nReading circuit declaration...")
                id_word += 1
                id_word = self.read_circuit_declaration(data_words, id_word)
                if id_word is None:
                    raise Exception("End of the circuit declaration section is not found.")
                continue

            if mix.compare_two_strings(word, "CIRCUIT_STRUCTURE"):
                id_word += 1
                id_word = self.read_circuit_structure(data_words, id_word)
                if id_word is None:
                    raise Exception("End of the circuit structure section is not found.")
                continue

            if mix.compare_two_strings(word, "LAUNCH_SECTION"):
                id_word += 1
                id_word = self.read_launch(data_words, id_word)

            # otherwise read next word:
            id_word += 1

    # -----------------------------------------------------------
    # --- REMOVE COMMENTS AND EMPTY LINES FROM AN INPUT FILE ---
    # -----------------------------------------------------------
    def remove_comments_empty_lines(self, data):
        data_lines = data.splitlines()
        data_lines_clr = []
        for one_line in data_lines:
            one_line = one_line.split(self.keyword_comment_, 1)[0]
            one_line = one_line.strip() # ro remove spaces
            data_lines_clr.append(one_line.lower())

        data_clean = mix.create_text_from_list(data_lines_clr)
        self.print_log("--- Original INPUT FILE ---", 0, True)
        self.print_log(data_clean, 1, True)
        self.print_log("----------------------------------", 0, True)

        data_words = []
        for one_line in data_lines_clr:
            for one_word in one_line.split():
                data_words.append(one_word)

        return data_words

    # -----------------------------------------------------------
    # --- READ CONSTANTS ---
    # -----------------------------------------------------------
    def read_constants(self, data_list, id_word):
        while id_word < len(data_list):
            word = data_list[id_word]
            if mix.compare_two_strings(word, "END_CONSTANTS"):
                id_word += 1
                return id_word

            # name of a constant:
            constant_name = word

            # read a value of the constant:
            id_word += 1
            word = data_list[id_word]
            constant_value = self.get_float(
                word, 
                f"A wrong format of the constant {constant_name} is given."
            )

            # save the constant:
            self.constants_[constant_name] = constant_value

            # write down the constant
            self.print_log(f"{constant_name} = {constant_value}")

            # read next word
            id_word += 1
        return None

    # -----------------------------------------------------------
    # --- GET AN INTEGER ---
    # -----------------------------------------------------------
    def get_integer(self, word, line_exc):
        res = re.search("<(.+?)>", word)
        if res is not None:
            constant_name = res.group(1)
            if constant_name in self.constants_:
                res_int = int(self.constants_[constant_name])
            else:
                raise Exception(f"A constant {constant_name} is not found.")
        else:
            if not mix.is_digit(word):
                raise Exception(line_exc)
            res_int = int(word)
        return res_int

    # -----------------------------------------------------------
    # --- GET A FLOAT ---
    # -----------------------------------------------------------
    def get_float(self, word, line_exc):
        res = re.search("<(.+?)>", word)
        if res is not None:
            constant_name = res.group(1)
            if constant_name in self.constants_:
                res_float = float(self.constants_[constant_name])
            else:
                raise Exception(f"A constant {constant_name} is not found.")
        else:
            try:
                res_float = float(word)
            except:
                raise Exception(line_exc)
        return res_float

    # -----------------------------------------------------------
    # --- READ DECLARATION OF CIRCUITS ---
    # -----------------------------------------------------------
    def read_circuit_declaration(self, data_list, id_word):
        while id_word < len(data_list):
            word = data_list[id_word]
            if mix.compare_two_strings(word, "END_CIRCUITS_DECLARATION"):
                id_word += 1
                return id_word

            # circuit name
            circ_name = word

            # save initial data of this circuit
            if circ_name in self.circuits_:
                self.print_log_warning(
                    f"A circuit {circ_name} has been already declared. " +
                    "The first declaration is taken."
                )
            else:
                self.circuits_[circ_name] = {
                    "flag_draw": False,
                    "flag_redraw": False,
                    "name": circ_name
                }

            # read registers of the circuit:
            id_word = self.read_reg_declaration(data_list, id_word, circ_name)

            # read next word
            id_word += 1
        return None

    # -----------------------------------------------------------
    # --- READ DESCRIPTION OF REGISTERS ---
    # -----------------------------------------------------------
    def read_reg_declaration(self, data_list, id_word, circ_name):
        dd_circ = self.circuits_[circ_name]

        # number of registers
        id_word += 1
        word = data_list[id_word]
        n_regs = self.get_integer(
            word, 
            f"Error: a wrong format of a number of registers for a circuit {circ_name}."
        )
        if n_regs < 1:
            raise Exception(f"Error: Circuit {circ_name}: every circuit must have at least one register.")

        # description of registers if present
        regs = {}
        qcr = None
        obj_reg = None
        names_regs = [None] * n_regs
        nq = 0 # total number of qubits in the circuit
        for id_reg in range(n_regs):

            # name of a register
            id_word += 1
            name_reg = data_list[id_word]
            names_regs[id_reg] = name_reg

            # number of qubits in a register
            id_word += 1
            word = data_list[id_word]
            nq_reg = self.get_integer(
                word, 
                f"Wrong format of a number of qubits in a register {name_reg}."
            )
            nq += nq_reg

            # ancilla flag (just skip it)
            id_word += 1
            word = data_list[id_word]

            # create a register
            obj_reg = qiskit.QuantumRegister(nq_reg, name_reg)

            # save the register:
            regs.update({
                name_reg: {
                    "name": name_reg,
                    "nq": nq_reg,
                    "obj": obj_reg,
                    "location": id_reg,
                }
            })

        # total number of qubits in the circuit:
        dd_circ["nq"] = nq

        # position of the first qubit within the register:
        shift_from_top = 0
        for id_reg in range(n_regs):
            name_reg = names_regs[id_reg]
            shift_from_top += regs[name_reg]["nq"]
            regs[name_reg]["id_first_qubit"] = nq - shift_from_top 

        # add registers to the quantum circuit:
        # !!! here, location of register is reversed  !!!
        # !!! to have correct visualisation in qiskit !!!
        for id_reg in reversed(range(n_regs)):
            obj_reg = regs[names_regs[id_reg]]["obj"]
            if qcr is None:
                qcr = qiskit.QuantumCircuit(obj_reg)
            else:
                qcr.add_register(obj_reg)

        # save the initial structure of the circuit:
        dd_circ.update({
            "n_regs": n_regs,
            "regs": regs,
            "obj": qcr
        })

        # output initial structure of the circuit: 
        self.print_log(f"...Circuit {circ_name}...")
        self.print_log("n of qubits: {:d}".format(dd_circ["nq"]), 1)
        self.print_log("n of registers: {:d}".format(dd_circ["n_regs"]), 1)     
        for name_reg in dd_circ["regs"].keys():
            nq = dd_circ["regs"][name_reg]["nq"]
            word_qu = "qubit" if nq == 1 else "qubits"
            self.print_log("register: {} with {:d} {}".format(name_reg, nq, word_qu), 2)

        return id_word

    # -----------------------------------------------------------
    # --- READ STRUCTURE OF A CIRCUIT ---
    # -----------------------------------------------------------
    def read_circuit_structure(self, data_list, id_word): 
        # circuit name
        circ_name = data_list[id_word]
        self.print_log(f"\nReading structure of a circuit {circ_name}...")
        if circ_name not in self.circuits_:
            self.print_log_warning(f"Circuit {circ_name} has not been declared. Skip it.")
            return id_word + 1
        dd_circ = self.circuits_[circ_name]
        qcr = dd_circ["obj"]

        # --- registers are organized from the high to low priority ---
        # (from top to bottom qubits):
        id_word += 1
        while id_word < len(data_list):
            word = data_list[id_word]

            if mix.compare_two_strings(word, "CIRCUIT_STRUCTURE"):
                return None

            if mix.compare_two_strings(word, "END_CIRCUIT_STRUCTURE"):
                id_word += 1
                return id_word

            if mix.compare_two_strings(word, "gate"):
                id_word += 1
                id_word = self.read_gate(data_list, id_word, qcr, dd_circ)
                continue

            if mix.compare_two_strings(word, "igate"):
                id_word += 1
                id_word = self.read_gate(data_list, id_word, qcr, dd_circ, True)
                continue

            if mix.compare_two_strings(word, "circuit"):
                id_word += 1
                id_word = self.read_subcircuit(data_list, id_word, qcr, dd_circ["nq"], circ_name)
                continue

            if mix.compare_two_strings(word, "icircuit"):
                id_word += 1
                id_word = self.read_subcircuit(data_list, id_word, qcr, dd_circ["nq"], circ_name, True)
                continue
        
            # otherwise read next word
            id_word += 1
        return None

    # -----------------------------------------------------------
    # --- READ a GATE ---
    # -----------------------------------------------------------
    def read_gate(self, data_list, id_word, qcr, dd_circ, flag_inv=False):
        def check_control(ids_target, ids_cs, ids_x, ids_c_it, ids_x_it):
            if ids_c_it:
                if not mix.is_list_rows_equal(ids_c_it):
                   raise Exception(
                        f"Circuit {circ_name}: Gate {gate_name}: " +
                        "the gate has several sets of iterative control nodes: they must be of the same size."
                    ) 
                if len(ids_target) is not len(ids_c_it[0]):
                    raise Exception(
                        f"Circuit {circ_name}: Gate {gate_name}: " +
                        "a number of nodes in one iterative control set must be equal to the number of target qubits."
                    )
            for id_target in ids_target:
                if id_target in ids_cs:
                    raise Exception(
                        f"Circuit {circ_name}: Gate {gate_name}: " +
                        "gate cannot have the same target and control qubits."
                    )
            return


        circ_name = dd_circ["name"]
        gate_name = data_list[id_word]

        # --- single qubit gates without parameteres ---
        if mix.is_string_among(gate_name, ["X", "Z", "H"]):
            # read a target qubit
            id_word, ids_target = self.read_target_qubit(data_list, id_word, gate_name, dd_circ)

            # read control qubits or end of the gate definition:
            id_word, flag_control, ids_cs, ids_x, ids_c_it, ids_x_it = \
                self.read_end_gate(data_list, id_word, gate_name, dd_circ)
            check_control(ids_target, ids_cs, ids_x, ids_c_it, ids_x_it)

            # add the gate to the circuit:
            count_id_target = -1
            for id_target in ids_target:
                count_id_target += 1
                if flag_control:
                    if mix.compare_two_strings(gate_name, "X"): oo_gate = qiskit.circuit.library.XGate()
                    if mix.compare_two_strings(gate_name, "Z"): oo_gate = qiskit.circuit.library.ZGate()
                    if mix.compare_two_strings(gate_name, "H"): oo_gate = qiskit.circuit.library.HGate()
                    self.set_single_qubit_gate(qcr, oo_gate, id_target, count_id_target, ids_cs, ids_x, ids_c_it, ids_x_it)
                else:
                    if mix.compare_two_strings(gate_name, "X"): qcr.x(id_target)
                    if mix.compare_two_strings(gate_name, "Z"): qcr.z(id_target)
                    if mix.compare_two_strings(gate_name, "H"): qcr.h(id_target)
            return id_word

        # --- single qubit gates with one parameter ---
        if mix.is_string_among(gate_name, ["Ry", "Rz", "Phase"]):
            # read a target qubit:
            id_word, ids_target = self.read_target_qubit(data_list, id_word, gate_name, dd_circ)

            # read the gate parameter:
            id_word, gate_parameter = self.read_gate_parameter(data_list, id_word, gate_name)

            # read control qubits or end of the gate definition:
            id_word, flag_control, ids_cs, ids_x, ids_c_it, ids_x_it = \
                self.read_end_gate(data_list, id_word, gate_name, dd_circ)

            count_id_target = -1
            for id_target in ids_target:
                count_id_target += 1
                if flag_control:
                    if mix.compare_two_strings(gate_name, "Ry"):    oo_gate = qiskit.circuit.library.RYGate(gate_parameter)
                    if mix.compare_two_strings(gate_name, "Rz"):    oo_gate = qiskit.circuit.library.RZGate(gate_parameter)
                    if mix.compare_two_strings(gate_name, "Phase"): oo_gate = qiskit.circuit.library.PhaseGate(gate_parameter)
                    self.set_single_qubit_gate(qcr, oo_gate, id_target, count_id_target, ids_cs, ids_x, ids_c_it, ids_x_it)
                else:
                    if mix.compare_two_strings(gate_name, "Ry"):    qcr.ry(gate_parameter, id_target)
                    if mix.compare_two_strings(gate_name, "Rz"):    qcr.rz(gate_parameter, id_target)
                    if mix.compare_two_strings(gate_name, "Phase"): qcr.p(gate_parameter, id_target)
            return id_word

        # --- Rc gate ---
        if mix.is_string_among(gate_name, ["Rc"]):
            # read a target qubit:
            id_word, ids_target = self.read_target_qubit(data_list, id_word, gate_name, dd_circ)

            # read the first parameter (angle for Rz):
            id_word, gate_parameter_1 = self.read_gate_parameter(data_list, id_word, gate_name)

            # read the second parameter (angle for Ry):
            id_word, gate_parameter_2 = self.read_gate_parameter(data_list, id_word, gate_name)

            # read control qubits or end the gate definition:
            id_word, flag_control, ids_cs, ids_x, ids_c_it, ids_x_it = \
                self.read_end_gate(data_list, id_word, gate_name, dd_circ)

            count_id_target = -1
            for id_target in ids_target:
                count_id_target += 1
                if flag_control:
                    oo_gate = qiskit.circuit.library.RZGate(gate_parameter_1)
                    oo_gate = qiskit.circuit.library.RYGate(gate_parameter_2)
                    self.set_single_qubit_gate(qcr, oo_gate, id_target, count_id_target, ids_cs, ids_x, ids_c_it, ids_x_it)
                else:
                    qcr.rz(gate_parameter_1, id_target)
                    qcr.ry(gate_parameter_2, id_target)
            return id_word

        # --- Conditional rotation ---
        if mix.is_string_among(gate_name, ["condR"]):
            # name of the conditional rotation:
            id_word += 1
            name_condr = data_list[id_word]
            if flag_inv: name_condr += "*"

            # read a target qubit:
            id_word, ids_target = self.read_target_qubit(data_list, id_word, gate_name, dd_circ)

            # read conditional nodes
            id_word, ids_conds = self.read_integer_for_qubit_positions(id_word, data_list, gate_name, dd_circ)

            # read control qubits or end of the gate definition:
            id_word, flag_control, ids_cs, ids_x, ids_c_it, ids_x_it = \
                self.read_end_gate(data_list, id_word, gate_name, dd_circ) 

            # create a box-gate to display the rotational gate:
            ids_qubits_basic = ids_conds + ids_cs
            ids_qubits_basic.sort()
            nq_gate = len(ids_qubits_basic) + 1
            for id_target in ids_target:
                ids_qubits = ids_qubits_basic + [id_target]
                ids_qubits.sort()

                qc_box = qiskit.QuantumCircuit(nq_gate, name=name_condr)
                gate_box = qc_box.to_gate()
                for id_x in ids_x: qcr.x(id_x)
                qcr.append(gate_box, ids_qubits)
                for id_x in ids_x: qcr.x(id_x)
            return id_word

        # --- Adder by one ---
        if mix.is_string_among(gate_name, ["adder1"]):
            # read target qubits:
            id_word, ids_target = self.read_target_qubit(data_list, id_word, gate_name, dd_circ)

            # read control qubits or end of the gate definition:
            id_word, flag_control, ids_cs, ids_x, ids_c_it, ids_x_it = \
                self.read_end_gate(data_list, id_word, gate_name, dd_circ) 

            # put the high-priority qubits at the beginning:
            ids_target.sort(reverse=True)
            nt = len(ids_target)

            # add CNOT and X gates with control nodes
            for id_x in ids_x: qcr.x(id_x)

            for id_t in range(nt-1):
                ids_cnot_cs = ids_target[id_t + 1:nt] + ids_cs
                qcr.mcx(ids_cnot_cs, ids_target[id_t])
            if ids_cs: qcr.mcx(ids_cs, ids_target[-1])
            else: qcr.x(ids_target[-1])

            for id_x in ids_x: qcr.x(id_x)

            return id_word

        # --- Subtractor by one ---
        if mix.is_string_among(gate_name, ["subtractor1"]):
            # read target qubits:
            id_word, ids_target = self.read_target_qubit(data_list, id_word, gate_name, dd_circ)

            # read control qubits or end of the gate definition:
            id_word, flag_control, ids_cs, ids_x, ids_c_it, ids_x_it = \
                self.read_end_gate(data_list, id_word, gate_name, dd_circ) 

            # put the low-priority qubits at the beginning:
            ids_target.sort(reverse=False)
            nt = len(ids_target)

            # add CNOT and X gates with control nodes
            for id_x in ids_x: qcr.x(id_x)

            if ids_cs: qcr.mcx(ids_cs, ids_target[0])
            else: qcr.x(ids_target[0])
            for id_t in range(1,nt):
                ids_cnot_cs = ids_target[0:id_t] + ids_cs
                qcr.mcx(ids_cnot_cs, ids_target[id_t])
            
            for id_x in ids_x: qcr.x(id_x)

            return id_word

        # --- SWAP ---
        if mix.is_string_among(gate_name, ["swap"]):
            # read target qubits:
            id_word, ids_target_1 = self.read_target_qubit(data_list, id_word, gate_name, dd_circ)
            id_word, ids_target_2 = self.read_target_qubit(data_list, id_word, gate_name, dd_circ)

            nt = len(ids_target_1)
            if np.abs(nt - len(ids_target_2)) > 0:
                raise Exception(
                            f"Circuit {circ_name}: Gate {gate_name}: " +
                            "number of source qubits must be equal to a number of target qubits."
                        )

            # read control qubits or end of the gate definition:
            id_word, flag_control, ids_cs, ids_x, ids_c_it, ids_x_it = \
                self.read_end_gate(data_list, id_word, gate_name, dd_circ) 

            # add CNOT gates
            for id_x in ids_x: qcr.x(id_x)
            for id_t in range(nt):
                t1 = ids_target_1[id_t]
                t2 = ids_target_2[id_t]
                qcr.mcx([t2] + ids_cs, t1)
                qcr.mcx([t1] + ids_cs, t2)
                qcr.mcx([t2] + ids_cs, t1)
            for id_x in ids_x: qcr.x(id_x)

            return id_word

        raise Exception(f"Gate {gate_name} is not found.")

    # -----------------------------------------------------------
    # --- Set a single-qubit gate ---
    # ----------------------------------------------------------
    def set_single_qubit_gate(self, qcr, oo_gate, id_target, count_id_target, ids_cs, ids_x, ids_c_it, ids_x_it):
        ids_c_it_for_one_target = []
        ids_x_it_for_one_target = []
        if ids_c_it:
            for ids_c_set in ids_c_it:
                ids_c_it_for_one_target += [ids_c_set[count_id_target]]
            for ids_x_set in ids_x_it:
                ids_x_it_for_one_target += [ids_x_set[count_id_target]]
        ids_cs_res = ids_cs + ids_c_it_for_one_target 
        ids_x_res = ids_x + ids_x_it_for_one_target

        ids_res_qubits = ids_cs_res + [id_target]

        for id_x in ids_x_res: qcr.x(id_x)
        qcr.append(oo_gate.control(len(ids_cs_res)), ids_res_qubits)
        for id_x in ids_x_res: qcr.x(id_x)
        return

    # -----------------------------------------------------------
    # --- READ a SUBCIRCUIT ---
    # -----------------------------------------------------------
    def read_subcircuit(self, data_list, id_word, qcr, nq, this_circ_name, flag_inv=False):
        name_subcircuit = data_list[id_word]
        print(f"Reading a subcircuit {name_subcircuit}")
        
        if name_subcircuit not in self.circuits_:
            self.print_log_warning(f"Circuit {name_subcircuit} is not found. Skip it.")
            return id_word
        if mix.compare_two_strings(name_subcircuit, this_circ_name):
            self.print_log_warning(f"Circuit {name_subcircuit} has itself as a subcircuit. Skip it")
            return id_word

        if flag_inv: name_subcircuit += "*"
        qc_box = qiskit.QuantumCircuit(nq, name=name_subcircuit)
        gate_box = qc_box.to_gate()

        all_qus = [i for i in range(nq)]
        qcr.append(gate_box, all_qus)
        return id_word

    # -----------------------------------------------------------
    # --- READ LAUNCH SECTION ---
    # -----------------------------------------------------------
    def read_launch(self, data_list, id_word):
        self.print_log("\nReading Launch Section...")
        while id_word < len(data_list):
            word = data_list[id_word]

            # end of the Launch Section
            if mix.compare_two_strings(word, "END_LAUNCH_SECTION"):
                return id_word

            # selector how to launch the program
            if mix.compare_two_strings(word, "sel_launch"):
                id_word += 1
                word = data_list[id_word]
                if not mix.is_string_among(word, self.available_launch_selectors_):
                    raise Exception(f"Selector {word} is unknown for the sel_launch.")
                self.sel_launch_ = word
                self.print_log(f"Launch selector: {self.sel_launch_}.")    
                id_word += 1
                continue

            # redraw a circuit after launching C++ QPlasma:        
            if mix.compare_two_strings(word, "flag_redraw"):
                id_word += 1
                name_circ = data_list[id_word]
                if name_circ in self.circuits_:
                    self.circuits_[name_circ]["flag_redraw"] = True
                    self.print_log(f"Circuit {name_circ} to redraw.")
                else:
                    self.print_log_warning(f"Circuit {name_circ} is not found. Skip it.")
                id_word += 1
                continue

            # draw a circuit:        
            if mix.compare_two_strings(word, "flag_draw"):
                id_word += 1
                name_circ = data_list[id_word]
                if name_circ in self.circuits_:
                    self.circuits_[name_circ]["flag_draw"] = True
                    self.print_log(f"Circuit {name_circ} to draw.")
                else:
                    self.print_log_warning(f"Circuit {name_circ} is not found. Skip it.")
                id_word += 1
                continue

            # unknown keyword
            id_word += 1
            self.print_log_warning(f"Unknown keyword {word} in the LAUNCH_SECTION. Skip it.")
        raise Exception("End of the LAUNCH_SECTION is not found.")

    # -----------------------------------------------------------
    # --- Read a target qubit ---
    # -----------------------------------------------------------
    def read_target_qubit(self, data_list, id_word, gate_name, dd_circ):
        circ_name = dd_circ["name"]
        nq = dd_circ["nq"]
        id_word, ids_target = self.read_integer_for_qubit_positions(id_word, data_list, gate_name, dd_circ, "target")
        for id_target in ids_target:
            if id_target >= nq:
                raise Exception(
                    f"In a definition of a gate {gate_name} of a circuit {circ_name}:\n" + 
                    f"the circuit has only {nq} qubits; id of a target qubit {id_target} is requested.\n"
                    f"id of a qubit can be only <({nq}-1) and >= 0."
                )
        return id_word, ids_target

    # ---------------------------------------------------------------
    # --- From register to qubit: find global location of a qubit ---
    # ---  via its local position within a register               ---
    # ---------------------------------------------------------------   
    def from_register_to_qubit(self, dd_reg, id_reg_qubit):
        id_qubit = dd_reg["id_first_qubit"] + id_reg_qubit
        return id_qubit

    # -----------------------------------------------------------
    # --- Read a parameter ---
    # -----------------------------------------------------------
    def read_gate_parameter(self, data_list, id_word, gate_name):
        id_word += 1
        word = data_list[id_word]
        gate_parameter = self.get_float(word, f"Wrong format of a gate parameter for a gate {gate_name}")
        return id_word, gate_parameter

    # -----------------------------------------------------------
    # --- END OF A GATE: read control nodes and               ---
    # ---                     end of the definition of a gate ---
    # -----------------------------------------------------------
    def read_end_gate(self, data_list, id_word, gate_name, dd_circ, key_gate="end_gate"):
        circ_name = dd_circ["name"]

        flag_control = False
        ids_c_res = []
        ids_x = []
        ids_c_it = []
        ids_x_it = []

        while id_word < len(data_list):
            id_word += 1
            word = data_list[id_word]

            if mix.compare_two_strings(word, key_gate):
                return id_word + 1, flag_control, ids_c_res, ids_x, ids_c_it, ids_x_it

            if mix.compare_two_strings(word, "control"):
                flag_control = True
                ids_c_res += self.read_control(data_list, id_word, gate_name, dd_circ)

            if mix.compare_two_strings(word, "ocontrol"):
                flag_control = True
                ids_x_one = self.read_control(data_list, id_word, gate_name, dd_circ)
                ids_c_res += ids_x_one
                ids_x     += ids_x_one

            if mix.compare_two_strings(word, "control_it"):
                flag_control = True
                ids_c_it_one = self.read_control(data_list, id_word, gate_name, dd_circ)
                ids_c_it += [ids_c_it_one]

            if mix.compare_two_strings(word, "ocontrol_it"):
                flag_control = True
                ids_x_it_one = self.read_control(data_list, id_word, gate_name, dd_circ)
                ids_c_it += [ids_x_it_one]
                ids_x_it += [ids_x_it_one]

        # wrong definition of a gate
        raise Exception(
            f"In a definition of a circuit {circ_name}: error in a definition of a gate {gate_name}."
        )

    # -----------------------------------------------------------
    # --- Read control nodes ---
    # -----------------------------------------------------------
    def read_control(self, data_list, id_word, gate_name, dd_circ):
        id_word, ids_c = self.read_integer_for_qubit_positions(id_word, data_list, gate_name, dd_circ)
        return ids_c

    # -----------------------------------------------------------
    # --- READ QUBIT POSITIONS AS AN INTEGER ---
    # -----------------------------------------------------------
    def read_integer_for_qubit_positions(self, id_word, data_list, gate_name, dd_circ, qu_type="control"):
        circ_name = dd_circ["name"]
        if dd_circ["n_regs"] == 0:
            raise Exception(
                f"A circuit {circ_name} does not have registers.\n" + 
                f"You cannot define the location of a gate {gate_name} via registers."
            )

        # number of registers:
        id_word += 1
        word = data_list[id_word]
        if not word.isdigit():
            n_regs = 1
            id_word -= 1
        else:
            n_regs = int(word)

        ids_c = []
        for id_reg in range(n_regs):
            id_word += 1
            reg_name = data_list[id_word]
            if reg_name not in dd_circ["regs"]:
                raise Exception(
                    f"In the definition of a gate {gate_name}:\n" + 
                    f"A circuit {circ_name} does not have a register with a name {reg_name}."
                )

            id_word += 1
            word = data_list[id_word]
            # if not mix.is_digit(word):
            #     raise Exception(
            #         f"In the definition of a gate {gate_name} in a circuit {circ_name}:\n"
            #         f"Wrong format of an integer for {qu_type} qubits in a register {reg_name}."
            #     )
            # integer_cs = int(word)

            line_exc = f"In the definition of a gate {gate_name} in a circuit {circ_name}:\n" + \
                f"Wrong format of an integer for {qu_type} qubits in a register {reg_name}."
            integer_cs = self.get_integer(word, line_exc)
            if integer_cs == 0:
                raise Exception(
                    f"In the definition of {qu_type} qubits in a register {reg_name} " + 
                    f"of a gate {gate_name} in a circuit {circ_name}:\n"
                    f"Integer cannot be =0"
                )

            nq_reg = dd_circ["regs"][reg_name]["nq"]
            max_integer = 2**nq_reg - 1
            if integer_cs < 0:
                integer_cs = max_integer + integer_cs + 1

            if integer_cs > max_integer:
                raise Exception(
                    f"In the definition of a gate {gate_name} in a circuit {circ_name}:\n"
                    f"Provided integer ({integer_cs}) for {qu_type} qubits " + 
                    f"is too large to be encoded in a register {reg_name}," + 
                    f"\nwhich has only {nq_reg} qubits (maximum possible integer is {max_integer})."
                )
            bit_array = mix.find_bit_array_of_int(integer_cs, nq_reg)

            # set control nodes where the bit array has 1:
            n_bitarray = len(bit_array)
            for id_bit in range(n_bitarray):
                one_bit = bit_array[n_bitarray - id_bit - 1]
                if one_bit == 1:
                    id_control = self.from_register_to_qubit(dd_circ["regs"][reg_name], id_bit)
                    ids_c.append(id_control)

        return id_word, ids_c

    # -----------------------------------------------------------
    # --- DRAW CIRCUITS ---
    # -----------------------------------------------------------
    def draw_circuits(self):
        print("\n*** DRAW CIRCUITS ***")
        for name_circuit in self.circuits_.keys():
            if self.circuits_[name_circuit]["flag_draw"]:
                qcr = self.circuits_[name_circuit]["obj"]

                # draw the circuit to a screen
                print(f"\nDraw a circuit {name_circuit}:")
                self.draw_one_circuit(qcr)

                # draw the circuit to a file
                filename_fig = self.path_ + "/" + self.name_ + "_" + \
                    name_circuit + self.fig_circuit_extention_
                self.draw_one_circuit(qcr, filename_fig)
        return

    # -----------------------------------------------------------
    # --- DRAW a CIRCUIT ---
    # -----------------------------------------------------------
    def draw_one_circuit(self, qcr, filename=None):
        circuit_length = mix.CIRCUIT_LENGTH

        if filename is not None:
            qcr.draw(output='mpl', filename=filename, fold=circuit_length, reverse_bits=True)
            return
        else:
            print(qcr.draw(fold=circuit_length, reverse_bits=True))
            return

    # -----------------------------------------------------------
    # --- LAUNCH a CIRCUIT ---
    # -----------------------------------------------------------
    def launch(self):
        if mix.compare_two_strings(self.sel_launch_, "draw"):
            return

        flag_test = "0"
        flag_compute_output = "0"
        flag_compute_iterator_output = "0"

        if mix.compare_two_strings(self.sel_launch_, "init_state"):
            flag_compute_output = "0"
            flag_compute_iterator_output = "0"
        if mix.compare_two_strings(self.sel_launch_, "full"):
            flag_compute_output = "1"
            flag_compute_iterator_output = "0"
        if mix.compare_two_strings(self.sel_launch_, "full_iterator"):
            flag_compute_output = "0"
            flag_compute_iterator_output = "1"

        # command_list = [
        #     "mpirun", "-np", str(self.n_mpi_procs_), 
        #     self.path_c_oracletool_,    # c++ oracletool executable file
        #     self.name_,
        #     self.path_,          # path to input files
        #     flag_compute_output, # do not compute output states only initial states
        #     flag_test            # flag to output additional information
        # ]
        command_list = [ 
            self.path_c_oracletool_,    # c++ oracletool executable file
            self.name_,
            self.path_,          # path to input files
            flag_compute_output, # whether compute output states from an oracle
            flag_test,            # flag to output additional information
            flag_compute_iterator_output # whether compute output states from an iterator
        ]
        print("\n\n=====================================================")
        print("=== Launch the c++ oracletool ===")
        print("=====================================================")
        print("Command list is ", command_list)
        res = subprocess.run(command_list)
        print("=== END OF C++ ORACLETOOL ===")
        print("=====================================================")

        if(np.abs(res.returncode) > 0):
            raise Exception(f"C++ oracletool finished with an error")

        # redraw some circuits after the analysis in the c++ oracletool:
        self.redraw_circuits()
        return

    # -----------------------------------------------------------
    # --- REDRAW CIRCUITS AFTER C++ ORACLETOOL ---
    # -----------------------------------------------------------
    def redraw_circuits(self):
        flag_first = True
        for name_circuit in self.circuits_.keys():
            flag_redraw = self.circuits_[name_circuit]["flag_redraw"]
            if not flag_redraw:
                continue

            if flag_first:
                flag_first = False
                print("\n*** REDRAW CIRCUITS ***")

            filename_fig = self.path_ + "/" + self.name_ + "_" + \
                name_circuit + "_redraw" + self.fig_circuit_extention_

            dc = dcm.DC()
            dc.path = self.path_
            dc.circuitFileName = name_circuit + self.circuit_ext_
            dc.filename_fig = filename_fig
            dc.draw()
        return

    # -----------------------------------------------------------
    # --- PRINT OUTPUT LOG DATA ---
    # -----------------------------------------------------------
    def print_log(self, data, int_indents=0, flag_only_file=False):
        data_write = data
        if int_indents > 0:
            line_indent = mix.LOG_INDENT * int_indents
            data_write = mix.insert_indent(
                data_write, line_indent
            )

        ff = open(self.log_name_, "a")
        ff.write(data_write + "\n")
        ff.close()
        if not flag_only_file:
            print(data_write)

    # -----------------------------------------------------------
    # --- PRINT OUTPUT LOG DATA WITH A WARNING ---
    # -----------------------------------------------------------
    def print_log_warning(self, data, int_indents=0):
        self.print_log("")
        self.print_log(f"WARNING: {data}\n", int_indents+1)


# -----------------------------------------------------------
# --- MAIN FUNCTION ---
# -----------------------------------------------------------
if __name__ == '__main__':
    oracletool = OracleTool()
    try:
        oracletool.read_input_file()
        oracletool.draw_circuits()
        oracletool.launch()
    except:
        e = sys.exc_info()
        oracletool.print_log(f"\nERROR: {e[1]}", 1)
        oracletool.close_oracletool()
        sys.exit(-1)
    
    oracletool.close_oracletool()
    sys.exit(0)







    




