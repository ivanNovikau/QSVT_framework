import qiskit
import math
import numpy as np

def reload():
    return

S_simulator = qiskit.Aer.backends(name='statevector_simulator')[0]

def analysis_prob(pp, qq):
    
    # good state:
    def is_good(q1):
        if q1[-1] == 1:
            return True
        else:
            return False
    
    pr_bad = 0 # probability of a bad state
    pr_good = 0 # probability of a good state
    for ii, p1 in enumerate(pp):
        if is_good(qq[ii]):
            pr_good += np.abs(p1)**2
        else:
            pr_bad += np.abs(p1)**2
    theta_a = np.arcsin(np.sqrt(pr_good))  # must be between 0 and np.pi/2
    n_opt = np.pi/(4*theta_a)

    print("prob. of GS: {:0.3e}".format(pr_good))
    print("prob. of BS: {:0.3e}".format(pr_bad))
    print("BS + GS: {:0.3e}".format(pr_bad + pr_good))
    print("optimal number of amplification operators: ", n_opt)

    return pr_good, pr_bad


def Wavefunction_adv(obj , *args, **kwargs):
    # Converts a number to binary, right to left LSB 152 153 o
    def Binary(number,total): 
        qubits = int(math.log(total,2))
        N = number
        b_num = np.zeros(qubits)
        for i in np.arange(qubits):
            if( N/((2)**(qubits-i-1)) >= 1 ):
                b_num[i] = 1
                N = N - 2 ** (qubits-i-1)
        B = [] 
        for j in np.arange(len(b_num)):
            B.append(int(b_num[j]))
        return B

    if(type(obj) == qiskit.QuantumCircuit ):
        statevec = qiskit.execute( obj, S_simulator, shots=1 ).result().get_statevector()
    if(type(obj) == np.ndarray):
        statevec = obj
        
    statevec = np.asarray(statevec)    
        
    sys = False
    NL = False
    dec = 5
    if 'precision' in kwargs:
        dec = int( kwargs['precision'] )
    if 'column' in kwargs:
        NL = kwargs['column']
    if 'systems' in kwargs:
        systems = kwargs['systems']
        sys = True
        last_sys = int(len(systems)-1)
        show_systems = []
        for s_chk in np.arange(len(systems)):
            if( type(systems[s_chk]) != int ):
                raise Exception('systems must be an array of all integers')
        if 'show_systems' in kwargs:
            show_systems = kwargs['show_systems']
            if( len(systems)!= len(show_systems) ):
                raise Exception('systems and show_systems need to be arrays of equal length')
            for ls in np.arange(len(show_systems)):
                if((show_systems[ls] != True) and (show_systems[ls] != False)):
                    raise Exception('show_systems must be an array of Truth Values')
                if(show_systems[ls] == True):
                    last_sys = int(ls) 
        else:
            for ss in np.arange(len(systems)):
                show_systems.append(True)

    if 'width' in kwargs:
        ww = kwargs["width"]  
        str_format = "{:<" + str(ww) + "s}"  

    if "ff" in kwargs:
        ff =  kwargs["ff"]     
    else:
        ff =  "{:20.3f}"

    str_wv = ''
    qubits = int(math.log(len(statevec),2))
    qq = []
    phases = []
    
    for i in range(len(statevec)):
        value = round(statevec[i].real, dec) + round(statevec[i].imag, dec) * 1j
        if np.abs(value) > 0:
            phases.append(value)
            
        if( (value.real != 0) or (value.imag != 0)):
            state = list(Binary(int(i),int(2**qubits)))
            state.reverse()
            state_str = ''
            qq.append(list(state))
            
            if( sys == True ): #Systems and SharSystems 
                k = 0 
                for s in np.arange(len(systems)):
                    if(show_systems[s] == True):
                        if(int(s) != last_sys):
                            state.insert(int(k + systems[s]), '>|' )
                            k = int(k + systems[s] + 1)
                        else:
                            k = int(k + systems[s]) 
                    else:
                        for s2 in np.arange(systems[s]):
                            del state[int(k)]
            for j in np.arange(len(state)):
                if(type(state[j])!= str):
                    state_str = state_str + str(int(state[j])) 
                else:
                    state_str = state_str + state[j]
            str_wv += ff.format(value) + ' |' + state_str + '>   '
            if(NL):
                str_wv = str_wv + '\n'

    return str_wv, phases, qq


def get_Fourier(nq, flag_also_inv=False, qc_fourier = None):
    m  = qiskit.QuantumRegister(nq, "m")

    if qc_fourier is None:
        qc_fourier = qiskit.QuantumCircuit(m, name="F")    
    for jj in range(nq):
        t = nq - 1 - jj
        qc_fourier.h(t)
        for kk in range(2, nq+1 - jj):
            oo_gate = qiskit.circuit.library.PhaseGate(2*np.pi/2**kk)
            oo_gate = oo_gate.control(1)
            all_qus = [t - (kk - 1)] + [t]
            qc_fourier.append(oo_gate, all_qus) 
    for jj in range(int(nq/2)):
        qc_fourier.swap(jj, nq-jj-1)

    gate_fourier = qc_fourier.to_gate()
    inv_gate_fourier = None
    if(flag_also_inv):
        inv_gate_fourier = gate_fourier.inverse()
    return qc_fourier, gate_fourier, inv_gate_fourier


def get_phase_estimation_circuit(nm, reg_target, gate_init, gate_u):
    m  = qiskit.QuantumRegister(nm, "m")
    cl = qiskit.ClassicalRegister(nm, "cl")
    qc_phase_estimation = qiskit.QuantumCircuit(m, reg_target, cl, name="PE")

    # prepare the initial state:
    qc_phase_estimation.append(gate_init, reg_target)

    # --- Set the initial superposition ---
    for ii in range(nm):
        qc_phase_estimation.h(m[ii])

    # --- add controlled U ---
    for ii in range(nm):
        nu = 2**(ii)
        qc_u = qiskit.QuantumCircuit(reg_target, name="U{:d}".format(nu))
        for _ in range(nu):
            qc_u.append(gate_u, reg_target)
        gate_uc = qc_u.to_gate().control(1)
        qc_phase_estimation.append(gate_uc, [m[ii]] + list(reg_target))
            
    # --- Inverse Fourier transform ---
    _, _, inv_gate_fourier = get_Fourier(nm, flag_also_inv=True)
    qc_phase_estimation.append(inv_gate_fourier, m)

    return qc_phase_estimation, m, cl


def get_phase_estimation(nm, reg_target, gate_init, gate_u):
    # res_phase = None
    # res_phase_next = None
    qc_phase_estimation = None

    qc_phase_estimation, m, cl = get_phase_estimation_circuit(nm, reg_target, gate_init, gate_u)

    # --- set the measurements --- 
    qc_phase_estimation.measure(m, cl)

    # --- perform the measurements ---
    backend = qiskit.BasicAer.get_backend('qasm_simulator')
    job = backend.run(qiskit.transpile(qc_phase_estimation, backend))
    data = job.result().get_counts(qc_phase_estimation)

    # --- find the full probability ---
    full_prob = 0
    for cl_state in data:   
        full_prob += data[cl_state]

    # --- find the estimated phase ---
    max_prob, max_prob_state, full_prob = 0, "", 0
    for cl_state in data:   
        prob_cl = data[cl_state]
        full_prob += prob_cl
        if max_prob < prob_cl:
            max_prob = prob_cl
            max_prob_state = cl_state
    max_prob_int = max_prob
    max_prob = (1. * max_prob) / full_prob
    res_int = string_bit_array_to_int(max_prob_state)
    res_phase = 2*np.pi * res_int / 2**nm

    next_prob, next_prob_state = 0, ""
    for cl_state in data:   
        prob_cl = data[cl_state]
        if next_prob < prob_cl and prob_cl < max_prob_int:
            next_prob = prob_cl
            next_prob_state = cl_state
    next_prob = (1. * next_prob) / full_prob
    res_int = string_bit_array_to_int(next_prob_state)
    res_phase_next = 2*np.pi * res_int / 2**nm

    # --- data to build histogram ---
    n_probs = len(data)
    hist_probs = np.zeros(n_probs)
    hist_phases = np.zeros(n_probs)
    count_hist = -1
    for cl_state in data:
        count_hist += 1
        int_1 = string_bit_array_to_int(cl_state)
        hist_phases[count_hist] = 2*np.pi * int_1 / 2**nm
        hist_probs[count_hist]  = (1.*data[cl_state]) / full_prob

    arg_sort = hist_phases.argsort()
    hist_phases = hist_phases[arg_sort]
    hist_probs = hist_probs[arg_sort]

    # --- results ---
    dd  = {
        "data": data,
        "circuit": qc_phase_estimation,
        "phase (max prob)": res_phase,
        "max. prob": max_prob,
        "phase (next prob)": res_phase_next,
        "next prob": next_prob,
        "max + next probs": max_prob + next_prob,
        "absolute error": np.abs(res_phase - res_phase_next) / 2.,
        "probs": hist_probs,
        "phases": hist_phases
    }
    return dd


def string_bit_array_to_int(string_line):
    res_int = 0
    nbits = len(string_line)
    counter_bit = nbits
    for one_char in string_line:
        counter_bit = counter_bit - 1
        if(int(one_char) == 1):
            res_int += 2**(counter_bit)
    return res_int
