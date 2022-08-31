#ifndef QCIRCUIT_H
#define QCIRCUIT_H

#include "QGates.h"
#include "CircuitLayers.h"

/**
 * @brief Circuit class.
 * ---
 * The 0-th qubit is the least significant 
 * (the rightmost in the statevector and the very bottom qubit in the quantum circuit)
 */
class QCircuit{
    public:
    /** Constructor of a circuit
     * @param[in] name name of the circuit;
     * @param[in] env input QuEST execution environment; 
     * @param[in] path_to_output path where output files should be written;
     * @param[in] nq number of qubits to create. 
     * If \p nq <= 0, the QuEST circuit will not be created. 
     * In this case, it will be necessary to launch separately the function "create";
     * @param[in] constants dictionary of constans used to create the circuit (by default, empty);
     * @param[in] flag_circuit print or not .circuit files (by default, false);
     * @param[in] flag_tex print or not .tex files (by default, false);
     * @param[in] flag_layers to calculate or not the layer for each gate;
     * */
    QCircuit(
        YCS name, const QuESTEnv& env, YCS path_to_output="./", YCU nq = 0,
        const std::map<std::string, qreal>& constants = std::map<std::string, qreal>(),
        YCB flag_circuit = false,
        YCB flag_tex = false,
        YCB flag_layers = false
    );

    /**
     * @brief Copy a circuit.
     * @param[in] oc is a circuit to copy.
     * @param[in] cname new name of the circuit.
     */
    QCircuit(YCCQ oc, YCS cname="");

    /** Destructor */
    ~QCircuit();

    /**
     * @brief Create a QuEST circuit.
     * @param nq number of qubits to create (has to be > 0);
     */
    void create(YCU nq);

    void create_circ_file();
    void create_tex_file();
    void finish_tex_file();

    void print_gates();

    inline std::string get_name(){ return name_; }
    inline unsigned get_n_qubits() const { return nq_; }
    inline int get_n_gates(){ return gates_.size(); }
    std::map<std::string, YVIv> get_regs(){ return regs_; }
    YVSv get_reg_names() const { return regnames_; }
    unsigned get_nq_in_reg(YCS reg_name) const{ return regs_.at(reg_name).size(); }

    /**
     * @brief Generate the circuit.
     */
    void generate();

    /**
     * @brief Generate the circuit taking into account stop gates.
     */
    void generate(std::string& stop_name, int& id_current);

    /**
     * @brief Get a conjugate transpose matrix.
     */
    void conjugate_transpose();

    inline Qureg get_qureg(){return c_;}

    /**
     * @brief Copy gates from the circuit \p circ, correct their positions \p regs_new, 
     * and add them to the end of the current circuit.
     * @param circ circuit from where new gates are taken;
     * @param regs_new registers of the current circuit to which the new gates are applied.
     *  \p regs_new.size() has to be = \p circ.nq_ .
     *  E.g., \p circ has two qubits: 0, 1 (0 - least significant (bottom) qubit).
     *  If \p regs_new = [3,2], then 0->3, 1->2: zero qubit of \p circ is applied to 
     *  the qubit 3 of the current circuit, while qubit 1 of \p circ is applied to
     *  the qubit 2 of the current circuit.
     * @param box a ghost gate to indicate boundaries where 
     *  gates from the circuit \p circ are placed into the current circuit.
     * @param flag_inv if true, first get conjugate transpose gates from \p circ. 
     * @param cs control qubits that should control each gate from \p circ.
     */
    void copy_gates_from(YCCQ circ, YCVI regs_new, 
        YCCB box = std::shared_ptr<const Box__>(nullptr),
        YCB flag_inv = false,
        YCVI cs = YVIv {}
    );

    /**
     * @brief Insert gates (!!! without copying them and without the correction of their qubits !!!) 
     * from the circuit \p of to the end of the current circuit.
     * @param of circuit from where new gates are taken;
     */
    void insert_gates_from(const QCircuit* of, 
        YCCB box = std::shared_ptr<const Box__>(nullptr));

    /**
     * @brief Add a register. 
     * The first added register is placed at the top of the circuit 
     * (is the most significant register).
     * Within the register, the 0-th elementin
     * @param[in] name name of the register;
     * @param[in] n_qubits number of qubits in the register;
     * @param[in] flag_ancilla if true, it is an ancilla register;
     * @return qubits of the register. 
     * The 0-th element in the returned vector corresponds to the least significant qubit.
     */
    YVIv add_register(YCS name, YCU n_qubits, YCB flag_ancilla=false);

    void set_standart_output_format();

    inline void set_ancillae(YCVI ancs){ ancs_ = YVIv(ancs); }

    /**
     * @brief Get the ancillae qubits of the circuit.
     * The vector with ancillae is not ordered.
     * @return YVIv unordered vector with ancillae qubits of the circuit.
     */
    inline YVIv get_ancillae() const { return ancs_; }

    /**
     * @brief Get an output format that corresponds to 
     *        the sizes of implied registers of the circut.
     * @return YVUv - vector that represents the format.
     */
    inline YVUv get_standart_output_format(){ return standart_output_format_; }

    /**
     * @brief Print locations of the qubit registers.
     */
    void print_reg_positions(std::ofstream& of) const;

    /**
     * @brief Save the qubit registers into the .circuit and .tex files.
     */
    void save_regs();

    /** @return pointer to the circuit. */
    YQCP get_the_circuit();

    /**
     * @brief Set the initial binary state of the circuit.
     */
    void set_init_binary_state(const bool& flag_mpi_bcast = false);

    void reset_init_vector(INIT_STATE__& state);

    /**
     * @brief Set initial amplitudes to specified qubits.
     * Elements \p ampl_vec_ in the state vector are filled starting from 
     * the elements corresponding to the low-priority qubit.
     * \p ampl_vec_real and \p ampl_vec_imag are assumed of the same size
     * @param[in] ampl_vec_real vector with real parts of amplitudes of size 2^nq;
     * @param[in] ampl_vec_imag vector with imaginary parts of amplitudes of size 2^nq;
     */
    void set_init_vector(YVQ ampl_vec_real, YVQ ampl_vec_imag);

    /**
     * @brief Set the qubit with id \p id_q to 1.
     * @param id_q id of the qubit to set to 1.
     */
    void set_qubit_state(YCU id_q);
    void set_qubit_state(YCVI ids_qs);

    /**
     * @brief Reset initial state of the circuit and remove all gates..
     */
    void reset();

    /**
     * @brief Reset only Qureg.
     */
    void reset_qureg();

    /**
     * @brief Empty saved binary initial states.
     */
    void empty_binary_states();

    /**
     * @brief Set zero state as an initial state of the circuit.
     */
    void prepare_zero_init_state();

    /**
     * @brief Set to 1 one of the qubits of a register. 
     * @param[in] name name of the register;
     * @param[in] id_reg_qubit qubit of the register to set to 1.
     */
    void set_reg_state(YCS name, YCI id_reg_qubit);

    /**
     * @brief Set to 1 several qubits of a register. 
     * @param[in] name name of the register;
     * @param ids_reg_qubits qubits of the register to set to 1.
     */
    void set_reg_state(YCS name, YCVI ids_reg_qubits);

    void read_structure_gate(
        YISS istr, YVI ids_target, qreal& par_gate, 
        YVI ids_control, YVI ids_x, 
        YVVI ids_control_it, YVVI ids_x_it
    );
    void read_structure_gate(
        YISS istr, YVI ids_target, qreal& par_gate1, qreal& par_gate2,
        YVI ids_control, YVI ids_x, 
        YVVI ids_control_it, YVVI ids_x_it
    );
    
    void read_structure_gate_condR_split(YISS istr, YCS path_in, YCB flag_inv=false);
    void read_structure_gate_adder1(YISS istr, YCS path_in, YCB flag_inv=false);
    void read_structure_gate_subtractor1(YISS istr, YCS path_in, YCB flag_inv=false);
    void read_structure_gate_adder(YISS istr, YCS path_in, YCB flag_inv=false);
    void read_structure_gate_swap(YISS istr, YCS path_in, YCB flag_inv=false);
    void read_structure_gate_fourier(YISS istr, YCS path_in, YCB flag_inv=false);
    void read_structure_gate_phase_estimation(
        YISS istr, YCS path_in, std::map<std::string, YSQ>& ocs, YCB flag_inv
    );
    void read_structure_gate_qsvt(
        YISS istr, YCS path_in, std::map<std::string, YSQ>& ocs, YCB flag_inv, QSVT_pars& data
    );

    /**
     * @brief store qubit indices to the output vector \p ids_target.
     * The qubits start from the less singificant one.
     */ 
    void read_reg_int(YISS istr, YVI ids_target, YCS word_start=std::string());
    void read_end_gate(YISS istr, YVI ids_control, YVI ids_x, YVVI ids_control_it, YVVI ids_x_it);

    template<class TGate>
    inline bool read_structure(YCS gate_name, YISS istr, YCB flag_inv=false)
    {
        std::vector<int> ids_target, ids_control, ids_x;
        YVVIv ids_control_it, ids_x_it;
        if(YMIX::compare_strings(gate_name, TGate::name_shared_, std::vector<std::string> {"X", "Y", "Z", "H"}))
        {
            qreal par_gate = nan("1");
            read_structure_gate(istr, ids_target, par_gate, ids_control, ids_x, ids_control_it, ids_x_it);

            if(ids_control_it.empty())
            {
                x(ids_x);
                for(auto const& id_target: ids_target) add_sqg<TGate>(id_target, ids_control, flag_inv);
                x(ids_x);
            }
            // --- iterative control nodes ---
            else
            {
                x(ids_x);
                int count_target = -1;
                
                for(auto const& id_target: ids_target)
                {
                    ++count_target;

                    YVIv ids_x_it_for_one_target, ids_control_res;
                    create_control_nodes(
                        count_target, ids_control_it, ids_x_it, ids_control, 
                        ids_control_res, ids_x_it_for_one_target
                    );

                    // add X gates for ocontrol and add a single-target gate
                    x(ids_x_it_for_one_target);
                    add_sqg<TGate>(id_target, ids_control_res, flag_inv);
                    x(ids_x_it_for_one_target);
                } 
                x(ids_x);
            }
            return true;
        }
        return false;
    }

    template<class TGate>
    inline bool read_structure(YCS gate_name, YISS istr, qreal& par_gate, YCB flag_inv=false)
    {
        std::vector<int> ids_target, ids_control, ids_x;
        YVVIv ids_control_it, ids_x_it;
        if(YMIX::compare_strings(gate_name, TGate::name_shared_, std::vector<std::string> {"Rx", "Ry", "Rz", "Phase"}))
        {
            read_structure_gate(istr, ids_target, par_gate, ids_control, ids_x, ids_control_it, ids_x_it);
            if(ids_control_it.empty())
            {
                x(ids_x);
                for(auto const& id_target: ids_target) add_sq_rg<TGate>(id_target, par_gate, ids_control, flag_inv);
                x(ids_x);
            }
            // --- iterative control nodes ---
            else
            {
                x(ids_x);
                int count_target = -1;
                
                for(auto const& id_target: ids_target)
                {
                    ++count_target;

                    YVIv ids_x_it_for_one_target, ids_control_res;
                    create_control_nodes(
                        count_target, ids_control_it, ids_x_it, ids_control, 
                        ids_control_res, ids_x_it_for_one_target
                    );

                    // add X gates for ocontrol and add a single-qubit gate with a parameter
                    x(ids_x_it_for_one_target);
                    add_sq_rg<TGate>(id_target, par_gate, ids_control_res, flag_inv);
                    x(ids_x_it_for_one_target);
                } 
                x(ids_x);
            }
            
            return true;
        }
        return false;
    }

    template<class TGate>
    inline bool read_structure(YCS gate_name, YISS istr, qreal& par_gate1, qreal& par_gate2, YCB flag_inv=false)
    {
        std::vector<int> ids_target, ids_control, ids_x;
        YVVIv ids_control_it, ids_x_it;
        if(YMIX::compare_strings(gate_name, TGate::name_shared_, std::vector<std::string> {"Rc"}))
        {
            read_structure_gate(istr, ids_target, par_gate1, par_gate2, ids_control, ids_x, ids_control_it, ids_x_it);
            if(ids_control_it.empty())
            {
                x(ids_x);
                for(auto const& id_target: ids_target) add_sq_rg<TGate>(id_target, par_gate1, par_gate2, ids_control, flag_inv);
                x(ids_x);
            }
            // --- iterative control nodes ---
            else
            {
                x(ids_x);
                int count_target = -1;
                
                for(auto const& id_target: ids_target)
                {
                    ++count_target;

                    YVIv ids_x_it_for_one_target, ids_control_res;
                    create_control_nodes(
                        count_target, ids_control_it, ids_x_it, ids_control, 
                        ids_control_res, ids_x_it_for_one_target
                    );

                    // add X gates for ocontrol and add a single-qubit gate with a parameter
                    x(ids_x_it_for_one_target);
                    add_sq_rg<TGate>(id_target, par_gate1, par_gate2, ids_control_res, flag_inv);
                    x(ids_x_it_for_one_target);
                } 
                x(ids_x);
            }
            
            return true;
        }
        return false;
    }

    inline void create_control_nodes(
        YCI count_target, YCVVI ids_control_it, YCVVI ids_x_it, YCVI ids_control, 
        YVI ids_control_res, YVI ids_x_it_for_one_target
    ){
        // take a control node from each set: one set is one "control" or "ocontrol"
        YVIv ids_c_it_for_one_target;
        ids_x_it_for_one_target = YVIv {};
        for(auto const& ids_c_set: ids_control_it) 
            ids_c_it_for_one_target.push_back(ids_c_set[count_target]);
        for(auto const& ids_x_set: ids_x_it)
            ids_x_it_for_one_target.push_back(ids_x_set[count_target]);

        // insert iterative control nodes to other control nodes
        ids_control_res = YVIv(ids_control);
        ids_control_res.insert(
            ids_control_res.end(), 
            ids_c_it_for_one_target.begin(), 
            ids_c_it_for_one_target.end()
        );
    }

    // add a single-qubit gate with several control gates:
    template<class TGate>
    YQCP add_sqg(YCI t, YVIv cs = {}, YCB flag_inv=false)
    {
        YSG oo = std::make_shared<TGate>(t);
        if(flag_inv) oo->conjugate_transpose();

        if(!cs.empty())
            oo->add_control_qubits(cs);
        gates_.push_back(oo);

        if(flag_layers_) oo_layers_->add_gate(oo);

        return get_the_circuit();
    }

    // add a single-qubit gate with one parameter with several control nodes: 
    template<class TGate>
    YQCP add_sq_rg(YCI t, YCQR a, YVIv cs = {}, YCB flag_inv=false)
    {
        YSG oo = std::make_shared<TGate>(t, a);
        if(flag_inv) oo->conjugate_transpose();

        if(!cs.empty())
            oo->add_control_qubits(cs);
        gates_.push_back(oo);

        if(flag_layers_) oo_layers_->add_gate(oo);

        return get_the_circuit();
    }

    // add a single-qubit gate with two parameters and with several control nodes: 
    template<class TGate>
    YQCP add_sq_rg(YCI t, YCQR a1, YCQR a2, YVIv cs = {}, YCB flag_inv=false)
    {
        YSG oo = std::make_shared<TGate>(t, a1, a2);
        if(flag_inv) oo->conjugate_transpose();

        if(!cs.empty())
            oo->add_control_qubits(cs);
        gates_.push_back(oo);

        if(flag_layers_) oo_layers_->add_gate(oo);

        return get_the_circuit();
    }

    /** Set Pauli X gate at \p t target qubit. 
     * @param[in] t target qubit;
     * @param[in] cs control qubits;
     * @return pointer to the circuit.
     * */
    inline YQCP x(YCI t, YVIv cs = {}){ return add_sqg<X__>(t, cs); }
    YQCP x(YCVI ts, YVIv cs = {});

    /** Set Pauli Y gate at \p t target qubit. 
     * @param[in] t target qubit;
     * @param[in] cs control qubits;
     * @return pointer to the circuit.
     * */
    inline YQCP y(YCI t, YVIv cs = {}){ return add_sqg<Y__>(t, cs); }
    YQCP y(YCVI ts, YVIv cs = {});

    /** Set Hadamard gate at \p t target qubit.  
     * @param[in] t target qubit;
     * @param[in] cs control qubits;
     * @return pointer to the circuit. 
     */
    inline YQCP h(YCI t, YVIv cs = {}){ return add_sqg<H__>(t, cs); }
    YQCP h(YCVI ts, YVIv cs = {});

    /** Set Pauli Z gate at \p t target qubit. 
     * @param[in] t target qubit;
     * @param[in] cs control qubits;
     * @return pointer to the circuit.
     * */
    inline YQCP z(YCI t, YVIv cs = {}){ return add_sqg<Z__>(t, cs); }
    YQCP z(YCVI ts, YVIv cs = {});

    /** Set a Rx-rotation gate.   
     * @param[in] t target qubit;
     * @param[in] a angle to rotate on;
     * @return pointer to the circuit.
     * */
    inline YQCP rx(YCI t, YCQR a, YVIv cs = {}, YCB flag_inv = false){ return add_sq_rg<Rx__>(t, a, cs, flag_inv); }

    /** Set a Ry-rotation gate.   
     * @param[in] t target qubit;
     * @param[in] a angle to rotate on;
     * @return pointer to the circuit.
     * */
    inline YQCP ry(YCI t, YCQR a, YVIv cs = {}, YCB flag_inv = false){ return add_sq_rg<Ry__>(t, a, cs, flag_inv); }

    /** Set a Rz-rotation gate.   
     * @param[in] t target qubit;
     * @param[in] a angle to rotate on;
     * @return pointer to the circuit.
     * */
    inline YQCP rz(YCI t, YCQR a, YVIv cs = {}, YCB flag_inv = false){ return add_sq_rg<Rz__>(t, a, cs, flag_inv); }

    /** Set the Rc-rotation gate: Ry(ay).Rz(az)   
     * @param[in] t target qubit;
     * @param[in] az angle of the Rz-rotation;
     * @param[in] ay angle of the Ry-rotation;
     * @return pointer to the circuit.
     * */
    inline YQCP rc(YCI t, YCQR az, YCQR ay, YVIv cs = {}, YCB flag_inv = false){ return add_sq_rg<Rc__>(t, az, ay, cs, flag_inv); }

    /** Set a phase shift gate.   
     * @param[in] t target qubit;
     * @param[in] a shift angle;
     * @return pointer to the circuit.
     * */
    inline YQCP phase(YCI t, YCQR a, YVIv cs = {}, YCB flag_inv = false){ return add_sq_rg<Phase__>(t, a, cs, flag_inv); }

    /**
     * @brief Add a CNOT gate.
     * @param[in] c control qubit;
     * @param[in] t target qubit.
     */
    inline YQCP cnot(YCI c, YCI t){ return x(t, YVIv {c}); }

    /** @brief integer encoded to \p ts is incremented */
    YQCP adder_by_one(YCVI ts, YCVI cs, YCB flag_inv);

    /** @brief integer encoded to \p ts is decremented */
    YQCP subtractor_by_one(YCVI ts, YCVI cs, YCB flag_inv);

    /** @brief addition of two variables (v1 and v2) encoded to the registers \p ts1 and \p ts2;
     * the three registers must be of the same size;
     * the output sum (v1 + v2) is written to the qubits [ts2[:], ts3[0]].
     * The register ts3 must be initialized to the zero state.
     * The qubits ts3[1:] are used to store carry bits and are returned in the zero state.
     * @param flag_box if true, draw the operator as a box, not as a circuit;
     */
    YQCP adder(YCVI ts1, YCVI ts2, YCVI ts3, YCVI cs, YCB flag_inv = false, YCB flag_box = false);

    /** @brief Add a swap operator between qubits \p t1 and \p t2. */
    YQCP swap(YCI t1, YCI t2, YVIv cs = {});

    /** @brief Quantum Fourier circuit placed to the qubits \p ts and
     *         controlled by the qubits \p cs;
     * @param flag_box if true, draw the operator as a box, not as a circuit;
     */
    YQCP quantum_fourier(YCVI ts, YCVI cs, YCB flag_inv = false, YCB flag_box = false);


    /** @brief Phase estimation (PE) operator to find an eigenphase of the operator \p A, 
     *         which sits on the qubits \p ta.
     *         The eigenstate of the circuit is prepared by the operator \p INIT.
     *         The final estimation is written to the qubits \p ty.
     *         The whole PE operator can be controlled by the qubits \p cs. 
     * @param flag_box if true, draw the operator as a box, not as a circuit;
     */ 
    YQCP phase_estimation(
        YCVI ta, 
        const std::shared_ptr<const QCircuit>& A, 
        const std::shared_ptr<const QCircuit>& INIT,
        YCVI ty, 
        YCVI cs, 
        YCB flag_inv = false,
        YCB flag_box = false
    );


    /** @brief QSVT inversion of the matrix encoded by the oracle \p BE, which sits on qubits \p qs_be.
     * The QSVT single rotations are placed at the qubit \p a_qsvt[0].
     * The whole QSVT circuit is controlled by the qubits \p cs.
     * @param flag_box if true, draw the operator as a box, not as a circuit;
     */
    YQCP qsvt_matrix_inversion(
        QSVT_pars& data,
        YCVI a_qsvt,
        YCVI qs_be, 
        const std::shared_ptr<const QCircuit> BE,
        YCVI cs, 
        YCB flag_inv = false,
        YCB flag_box = false
    );


    /**
     * @brief Add a Stop gate to a quantum state at this point.
     * @param name is a name to identify this stop position.
     * @return YQCP 
     */
    inline YQCP add_stop_gate(YCS name)
    {
        YSG oo = std::make_shared<GStop__>(name);
        gates_.push_back(oo);
        return get_the_circuit();
    }

    inline void add_gate(std::shared_ptr<Gate__> oo){ gates_.push_back(oo); }

    /**
     * @brief Return the reference to the circuit statevector. 
     */
    void get_ref_to_state_vector(qreal*& state_real, qreal*& state_imag);

    /**
     * Return states with nonzero amplitudes.
     * @param[in] flag_ZeroHighPriorAnc if true, then assume that all ancillae are 
     * the high-priority qubits, and calculate only states where these ancillae are in the zero state.
     */
    void get_state(YMIX::StateVectorOut& out, YCB flag_ZeroHighPriorAnc = false);

    /**
     * @brief Make the entire circuit controlled by qubit cs.
     * @param cs new controlled qubits of the circuit.
     */
    void controlled(YCVI cs);
    void controlled(YCI c);

    inline std::string get_name() const { return name_; }

    /**
     * @brief Get a number of ancilla qubits in the circuit. 
     */
    inline uint32_t get_na() const {return ancs_.size();}

private:
    qreal get_value_from_word(YCS word);
    void  qsvt_read_parameters(YCS filename, QSVT_pars& data);



private:
    // if one addes a new property, do not forget to add it to the copy constructor.

    std::string name_; // name of the circuit;
    QuESTEnv env_; // QuEST execution environment;
    Qureg c_; // quantum register (circuit);
    YMIX::YTimer timer_; // timer for the circuit object;
    std::string path_to_output_; // path where output files to be written
    unsigned nq_; // number of qubits in the circuit;

    INIT_STATE__ init_state_; // initial state;

    std::string cfname_;  // name of the .circuit file;
    std::string texname_; // name of the .tex file;

    // registers:
    // regs_[rname][i] is the i-th qubit in the register "rname";
    // the 0-th qubit is the least significant in the register.
    std::map<std::string, YVIv> regs_; 

    // which registers are ancilla:
    std::map<std::string, bool> flags_anc_regs_;

    /** register names: the first name corresponds to the register at the top. */
    std::vector<std::string> regnames_; 

    // initial state of the circuit as an array of bits;
    // ib_state_[nq-1] corresponds to the least significant qubit 
    // (qubit at the bottom, the very right qubit)
    std::vector<short> ib_state_; 

    std::vector<YSG> gates_; // gates in the circuit;

    /** @brief Index of the starting gate to generate the circuit.*/
    int id_start_ = 0;

    std::vector<unsigned> standart_output_format_; // standart output format

    std::map<std::string, qreal> constants_;

    YVIv ancs_; // position of ancilla qubits (unordered vector);

    // object to organise gates in layers:
    std::shared_ptr<CircuitLayers__> oo_layers_;

    // circ_lines[i][j]: j-th phrase in the i-th line from the top:
    std::vector<std::vector<std::string>> tex_lines_;

    // tex_noc_[i] = id of the next non-occupied position on the i-th qubit:
    std::vector<uint64_t> tex_noc_;

    // print or not .circuit files:
    bool flag_circuit_;

    //print or not .tex files:
    bool flag_tex_;

    // to calculate or not the layers;
    bool flag_layers_; 
};






#endif