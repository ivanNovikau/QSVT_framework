import qiskit
from qiskit.tools.visualization import circuit_drawer
import sys
import numpy as np
import math

try:
    import pylib.mix as mix
except:
    import mix as mix


def reload():
    mix.reload_module(mix)
    return

class DC:
    path = None # path to a file, where a circuit of interest is described;
    circuitFileName = None # name of a file where a circuit of interest is described;
    gates = None  # list of maps with the description of gates 

    circuit = None  # quantum circuit;
    circuitName = None  # name of the circuit;
    nq = None  # number of qubits in the circuit

    qregs = None  # {name_reg: n_qubits}
    regs_nq = None
    qregs_names = None
    flags_reg_anc = None

    filename_fig = None

    circuit_fold = None

    def __init__(self, path_circ, name_circ):
        self.path = path_circ
        self.circuitFileName = name_circ + ".circuit"
        self.filename_fig = self.path + "/" + name_circ + ".png"
        self.gates = []  
        self.circuitName = ""  
        self.nq = None   
        self.qregs = {}  
        self.regs_nq = []
        self.qregs_names = []
        self.flags_reg_anc = []
        self.circuit_fold = mix.CIRCUIT_LENGTH


    def draw(self):
        print(self.circuit.draw(fold=self.circuit_fold, reverse_bits=True))
        self.circuit.draw(output='mpl', filename=self.filename_fig, fold=self.circuit_fold, reverse_bits=True)
        print(f"the circuit depth: {self.circuit.depth()}")
        return


    def add_mc_gate(self, oo_gate, gate, qcr):
        oo_gate = oo_gate.control(len(gate["controls"]))
        all_qus = gate["controls"] + [gate["targets"][0]]
        qcr.append(oo_gate, all_qus)


    # reverse the position of a qubit;
    # the QSVT_framework and qiskit have opposite ordering of qubits: 
    def get_q(self, id_q):
        # return (self.nq - id_q - 1)
        return id_q


    def read_default(self, counter_line, dg, gmap, flag_box = False):
        # --- name of the gate ---
        id_field = 0
        gmap["name"] = dg[id_field]
        if flag_box:
            gmap["name"] += "-box"

        # print("Reading gate " + gmap["name"])

        # --- whether the gate is complex-conjugated ---
        id_field += 1
        if mix.compare_two_strings(dg[id_field], "conj"):
            gmap["flag-conj"] = True 
        elif mix.compare_two_strings(dg[id_field], "orig"):
            gmap["flag-conj"] = False 
        else:
            raise Exception(
                "Line {:d}: after the gate name must be the flag ".format(counter_line) +
                "of the complex-conjugation: orig or conj."
            )

        # --- id of layer where the gate sits ---
        id_field += 1
        if not mix.compare_two_strings(dg[id_field], "layer"):
            raise Exception("Line {:d}: error in the layer id.".format(counter_line))
        id_field += 1
        gmap["id-layer"] = int(dg[id_field])

        # --- target qubits ---
        id_field += 1
        if not mix.compare_two_strings(dg[id_field], "targets"):
            raise Exception("Line {:d}: error in the target qubits.".format(counter_line))

        id_field += 1
        nt = int(dg[id_field])
        gmap["nt"] = nt

        id_field += 1
        gmap["targets"] = [0]*nt
        for i in range(nt):
            gmap["targets"][i] = self.get_q(int(dg[id_field+i]))

        # --- control qubits ---
        id_field += nt
        if not mix.compare_two_strings(dg[id_field], "controls"):
            raise Exception("Line {:d}: error in the control qubits.".format(counter_line))

        id_field += 1
        nc = int(dg[id_field])
        gmap["nc"] = nc

        id_field += 1
        gmap["controls"] = [0]*nc
        for i in range(nc):
            gmap["controls"][i] = self.get_q(int(dg[id_field+i]))

        # --- gate parameters ---
        id_field += nc
        if not mix.compare_two_strings(dg[id_field], "pars"):
            raise Exception("Line {:d}: error in the gate parameters.".format(counter_line))

        id_field += 1
        np = int(dg[id_field])
        gmap["np"] = np

        id_field += 1
        gmap["parameters"] = [0]*np
        for i in range(np):
            gmap["parameters"][i] = float(dg[id_field+i]) 

        # --- store the read gate ---
        self.gates.append(gmap)


    # Read the .circuit file:
    def read_file(self):
        import csv

        try:
            file_name = self.path + "/" + self.circuitFileName
            print("Reading the file " + file_name + "...\n")

            data = []
            with open(file_name, newline='') as csvfile:
                spamreader = csv.reader(csvfile, delimiter=' ', quotechar='|')
                for xx in spamreader:
                    data.append(xx)
                        
            # circuit name:
            self.circuitName = data[0][0]

            # number of qubits in the circuit:
            self.nq = int(data[0][1])

            # circuit registers:
            dg = data[1]
            if mix.compare_two_strings("QubitRegisterNames", dg[0]):
                for i in range(1,len(dg)-2,3):
                    self.qregs_names.append(dg[i])

                    nq_in_reg = int(dg[i+1])
                    self.regs_nq.append(nq_in_reg)
                    self.qregs[dg[i]] = nq_in_reg

                    self.flags_reg_anc.append(dg[i+2])
            else:
                raise Exception("Registers are not indicated in the .circuit file.")
                
            # read gates of the circuit:
            counter_line = 2
            gatesData = data[counter_line:len(data)]
            flag_box = False
            name_box = ""
            for dg in gatesData:
                gmap = {}
                gname = dg[0]
                counter_line += 1

                # initial states:
                if mix.compare_two_strings(gname, "initState"):
                    bb = np.zeros(self.nq)
                    for i in range(self.nq):
                        bb[i] = int(dg[i+1])
                    gmap["initState"] = bb
                    self.gates.append(gmap)
                    continue

                # do not read gates inside boxes:
                if flag_box and not mix.compare_two_strings(gname, "Box"):
                    continue
                
                # Box gate to combine a set of gates into one gate during the drawing:
                if mix.compare_two_strings(gname, "Box"):
                     # to take into account only the biggest boxes and 
                    #   do not include the interior ones:
                    if flag_box: 
                        if name_box != dg[2]:
                            continue

                    flag_box = True if mix.compare_two_strings(dg[1], "start") else False

                    # if it is the beginnig of the box, store only its name:
                    if flag_box:
                        name_box = dg[2]

                    # if it is the end of the box, store all its parameters:
                    if not flag_box:
                        self.read_default(counter_line, dg[2:], gmap, flag_box = True)
                    continue

                # other gates:
                self.read_default(counter_line, dg[0:], gmap, flag_box = False)

        except Exception as e: 
            print("--------------------------------")
            print(e)
            print("--------------------------------")
            

        print("Name of the circuit: ", self.circuitName)
        print("Number of qubits: ", self.nq)


    # create a QISKIT circuit:
    def form_circuit(self):
        qcr = None 
        regs_names_inv = np.flip(self.qregs_names)
        for reg_name in regs_names_inv:
            reg_size = self.qregs[reg_name]
            obj_reg = qiskit.QuantumRegister(reg_size, reg_name)
            if qcr is None:
                qcr = qiskit.QuantumCircuit(obj_reg)
            else:
                qcr.add_register(obj_reg)

        # add gates to the circuit:
        for gate in self.gates:
            name = gate["name"]

            # if mix.compare_two_strings(name, "initState"):
            #     bb = gate["initState"]
            #     for i in range(self.nq):
            #         if bb[i] == 1:
            #             qcr.initialize([0,1], i)

            if mix.compare_two_strings(name, "X"):
                if gate["nc"] > 0:
                    qcr.mcx(gate["controls"], gate["targets"][0])
                else:
                    qcr.x(gate["targets"][0])
                continue

            if mix.compare_two_strings(name, "Y"):
                if gate["nc"] > 0:
                    oo_gate = qiskit.circuit.library.YGate()
                    self.add_mc_gate(oo_gate, gate, qcr)
                else:
                    qcr.y(gate["targets"][0])
                continue

            if mix.compare_two_strings(name, "Z"):
                if gate["nc"] > 0:
                    oo_gate = qiskit.circuit.library.ZGate()
                    self.add_mc_gate(oo_gate, gate, qcr)
                else:
                    qcr.z(gate["targets"][0])
                continue

            if mix.compare_two_strings(name, "H"):
                if gate["nc"] > 0:
                    oo_gate = qiskit.circuit.library.HGate()
                    self.add_mc_gate(oo_gate, gate, qcr)
                else:
                    qcr.h(gate["targets"][0])
                continue

            if mix.compare_two_strings(name, "Rx"):
                par = gate["parameters"][0]
                if gate["flag-conj"]:
                    par = - par
                if gate["nc"] > 0:
                    oo_gate = qiskit.circuit.library.RXGate(par)
                    self.add_mc_gate(oo_gate, gate, qcr)
                else:
                    qcr.rx(par, gate["targets"][0])
                continue

            if mix.compare_two_strings(name, "Ry"):
                par = gate["parameters"][0]
                if gate["flag-conj"]:
                    par = - par
                if gate["nc"] > 0:
                    oo_gate = qiskit.circuit.library.RYGate(par)
                    self.add_mc_gate(oo_gate, gate, qcr)
                else:
                    qcr.ry(par, gate["targets"][0])
                continue

            if mix.compare_two_strings(name, "Rz"):
                par = gate["parameters"][0]
                if gate["flag-conj"]:
                    par = - par
                if gate["nc"] > 0:
                    oo_gate = qiskit.circuit.library.RZGate(par)
                    self.add_mc_gate(oo_gate, gate, qcr)
                else:
                    qcr.rz(par, gate["targets"][0])
                continue

            if mix.compare_two_strings(name, "Rc"):
                par_z = gate["parameters"][0]
                par_y = gate["parameters"][1]
                if gate["flag-conj"]:
                    par_z = - par_z
                    par_y = - par_y

                if gate["nc"] > 0:
                    oo_gate_z = qiskit.circuit.library.RZGate(par_z)
                    oo_gate_y = qiskit.circuit.library.RYGate(par_y)
                    if gate["flag-conj"]:
                        self.add_mc_gate(oo_gate_y, gate, qcr)
                        self.add_mc_gate(oo_gate_z, gate, qcr)
                    else:
                        self.add_mc_gate(oo_gate_z, gate, qcr)
                        self.add_mc_gate(oo_gate_y, gate, qcr)
                else:
                    if gate["flag-conj"]:
                        qcr.ry(par_y, gate["targets"][0])
                        qcr.rz(par_z, gate["targets"][0])
                    else:
                        qcr.rz(par_z, gate["targets"][0])
                        qcr.ry(par_y, gate["targets"][0])
                continue

            if mix.compare_two_strings(name, "Phase"):
                par = gate["parameters"][0]
                if gate["flag-conj"]:
                    par = - par
                if gate["nc"] > 0:
                    oo_gate = qiskit.circuit.library.PhaseGate(par)
                    self.add_mc_gate(oo_gate, gate, qcr)
                else:
                    qcr.p(par, gate["targets"][0])
                continue

            if "-box" in name:
                name_show = gate["name"][0:len(name_show)-4]
                qc_box = qiskit.QuantumCircuit(len(gate["targets"]), name=name_show)
                gate_box = qc_box.to_gate()
                if gate["nc"] > 0:
                    gate_box = gate_box.control(gate["nc"])

                temp = gate["targets"]
                temp.reverse()
                all_qus = gate["controls"] + temp
                qcr.append(gate_box, all_qus)
                continue

            print("Warning: do not know how to plot the gate with the name [{:s}]".format(name))

        # save the created circuit:
        self.circuit = qcr


    def get_bit_array(self, number): 
        N = number
        b_num = np.zeros(self.nq)
        for i in np.arange(self.nq):
            if( N/((2)**(self.nq-i-1)) >= 1 ):
                b_num[i] = 1
                N = N - 2 ** (self.nq-i-1)
        B = [] 
        for j in np.arange(len(b_num)):
            B.append(int(b_num[j]))
        return B


    def form_output_states(self, format_state = [], ff = "{:23.3e}"):
        if not format_state:
            format_state = [self.nq]

        # --- calculate the state vector ---
        statevec = qiskit.execute(
            self.circuit, 
            qiskit.Aer.backends(name='statevector_simulator')[0], 
            shots=1 
        ).result().get_statevector()   
        statevec = np.asarray(statevec)    
            
        # --- form the states and their amplitudes ---
        str_wv = ''
        qq = []
        ampls = []
        for i in range(len(statevec)):
            value = statevec[i].real + statevec[i].imag * 1j
            if mix.is_zero(value):
                continue
            else:
                ampls.append(value)
                
            if( (value.real != 0) or (value.imag != 0)):
                state = list(self.get_bit_array(int(i)))
                # state.reverse()
                state_str = ''
                qq.append(list(state))
                 
                k = 0 
                for s in np.arange(len(format_state)):
                    if(s != len(format_state) - 1):
                        state.insert(k + format_state[s], '>|' )
                        k = int(k + format_state[s] + 1)
                    else:
                        k = int(k + format_state[s]) 
                for j in np.arange(len(state)):
                    if(type(state[j])!= str):
                        state_str = state_str + str(int(state[j])) 
                    else:
                        state_str = state_str + state[j]
                str_wv += ff.format(value) + ' |' + state_str + '>   '
                str_wv = str_wv + '\n'
        return str_wv, ampls, qq






    # def draw(self):
    #     # read a file with a circuit description:
    #     self.read_file()

    #     # Quantum circuit:
    #     qcr = None 
    #     regs_names_inv = np.flip(self.qregs_names)
    #     for reg_name in regs_names_inv:
    #         reg_size = self.qregs[reg_name]
    #         obj_reg = qiskit.QuantumRegister(reg_size, reg_name)
    #         if qcr is None:
    #             qcr = qiskit.QuantumCircuit(obj_reg)
    #         else:
    #             qcr.add_register(obj_reg)

    #     # add gates to the circuit:
    #     for gate in self.gates:
    #         name = gate["name"]

    #         # if mix.compare_two_strings(name, "initState"):
    #         #     bb = gate["initState"]
    #         #     for i in range(self.nq):
    #         #         if bb[i] == 1:
    #         #             qcr.initialize([0,1], i)

    #         if mix.compare_two_strings(name, "X"):
    #             if gate["nc"] > 0:
    #                 qcr.mcx(gate["controls"], gate["targets"][0])
    #             else:
    #                 qcr.x(gate["targets"][0])
    #             continue

    #         if mix.compare_two_strings(name, "Y"):
    #             if gate["nc"] > 0:
    #                 oo_gate = qiskit.circuit.library.YGate()
    #                 self.add_mc_gate(oo_gate, gate, qcr)
    #             else:
    #                 qcr.y(gate["targets"][0])
    #             continue

    #         if mix.compare_two_strings(name, "Z"):
    #             if gate["nc"] > 0:
    #                 oo_gate = qiskit.circuit.library.ZGate()
    #                 self.add_mc_gate(oo_gate, gate, qcr)
    #             else:
    #                 qcr.z(gate["targets"][0])
    #             continue

    #         if mix.compare_two_strings(name, "H"):
    #             if gate["nc"] > 0:
    #                 oo_gate = qiskit.circuit.library.HGate()
    #                 self.add_mc_gate(oo_gate, gate, qcr)
    #             else:
    #                 qcr.h(gate["targets"][0])
    #             continue

    #         if mix.compare_two_strings(name, "Rx"):
    #             par = gate["parameters"][0]
    #             if gate["flag-conj"]:
    #                 par = - par
    #             if gate["nc"] > 0:
    #                 oo_gate = qiskit.circuit.library.RXGate(par)
    #                 self.add_mc_gate(oo_gate, gate, qcr)
    #             else:
    #                 qcr.rx(par, gate["targets"][0])
    #             continue

    #         if mix.compare_two_strings(name, "Ry"):
    #             par = gate["parameters"][0]
    #             if gate["flag-conj"]:
    #                 par = - par
    #             if gate["nc"] > 0:
    #                 oo_gate = qiskit.circuit.library.RYGate(par)
    #                 self.add_mc_gate(oo_gate, gate, qcr)
    #             else:
    #                 qcr.ry(par, gate["targets"][0])
    #             continue

    #         if mix.compare_two_strings(name, "Rz"):
    #             par = gate["parameters"][0]
    #             if gate["flag-conj"]:
    #                 par = - par
    #             if gate["nc"] > 0:
    #                 oo_gate = qiskit.circuit.library.RZGate(par)
    #                 self.add_mc_gate(oo_gate, gate, qcr)
    #             else:
    #                 qcr.rz(par, gate["targets"][0])
    #             continue

    #         if mix.compare_two_strings(name, "Rc"):
    #             par_z = gate["parameters"][0]
    #             par_y = gate["parameters"][1]
    #             if gate["flag-conj"]:
    #                 par_z = - par_z
    #                 par_y = - par_y

    #             if gate["nc"] > 0:
    #                 oo_gate_z = qiskit.circuit.library.RZGate(par_z)
    #                 oo_gate_y = qiskit.circuit.library.RYGate(par_y)
    #                 if gate["flag-conj"]:
    #                     self.add_mc_gate(oo_gate_y, gate, qcr)
    #                     self.add_mc_gate(oo_gate_z, gate, qcr)
    #                 else:
    #                     self.add_mc_gate(oo_gate_z, gate, qcr)
    #                     self.add_mc_gate(oo_gate_y, gate, qcr)
    #             else:
    #                 if gate["flag-conj"]:
    #                     qcr.ry(par_y, gate["targets"][0])
    #                     qcr.rz(par_z, gate["targets"][0])
    #                 else:
    #                     qcr.rz(par_z, gate["targets"][0])
    #                     qcr.ry(par_y, gate["targets"][0])
    #             continue

    #         if mix.compare_two_strings(name, "Phase"):
    #             par = gate["parameters"][0]
    #             if gate["flag-conj"]:
    #                 par = - par
    #             if gate["nc"] > 0:
    #                 oo_gate = qiskit.circuit.library.PhaseGate(par)
    #                 self.add_mc_gate(oo_gate, gate, qcr)
    #             else:
    #                 qcr.p(par, gate["targets"][0])
    #             continue

    #         if "-box" in name:
    #             name_show = gate["name"][0:len(name_show)-4]
    #             qc_box = qiskit.QuantumCircuit(len(gate["targets"]), name=name_show)
    #             gate_box = qc_box.to_gate()
    #             if gate["nc"] > 0:
    #                 gate_box = gate_box.control(gate["nc"])

    #             temp = gate["targets"]
    #             temp.reverse()
    #             all_qus = gate["controls"] + temp
    #             qcr.append(gate_box, all_qus)
    #             continue

    #         print("Warning: do not know how to plot the gate with the name [{:s}]".format(name))


    #     # draw the circuit:
    #     print(qcr.draw(fold=self.circuit_fold, reverse_bits=True))
    #     qcr.draw(output='mpl', filename=self.filename_fig, fold=self.circuit_fold, reverse_bits=True)
    #     # qcr.draw(filename=self.filename_fig, fold=self.circuit_fold, reverse_bits=True)

    #     # circuit depth:
    #     print(f"the circuit depth: {qcr.depth()}")
    #     return