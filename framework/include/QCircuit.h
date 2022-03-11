#ifndef QCIRCUIT_H
#define QCIRCUIT_H

#include "QGates.h"
#include "CircuitLayers.h"

/**
 * @brief Circuit class.
 * ---
 * The 0-th qubit is the least significant 
 * (the rightmost in the statevector and the very bottom qubit in the quantum circuit)
 * ---
 * Definitions:
 * -> CF is the text representation of a circuit component (e.g. gate) in the circuit file.
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
     * */
    QCircuit(
        YCS name, const QuESTEnv& env, YCS path_to_output="./", YCU nq = 0,
        const std::map<std::string, qreal>& constants = std::map<std::string, qreal>()
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

    void print_gates(YCB flag_print = false);

    inline std::string get_name(){ return name_; }
    inline unsigned get_n_qubits(){ return nq_; }
    inline int get_n_gates(){ return gates_.size(); }
    std::map<std::string, YVIv> get_regs(){ return regs_; }
    YVSv get_reg_names() const { return regnames_; }
    unsigned get_nq_in_reg(YCS reg_name) const{ return regs_.at(reg_name).size(); }

    /**
     * @brief Generate the circuit.
     */
    void generate(const bool& flag_print = false);

    /**
     * @brief Generate the circuit taking into account stop gates.
     */
    void generate(std::string& stop_name, int& id_current, const bool& flag_print = false);

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
     *  If \p circ = [3,2], then 0->3, 1->2: zero qubit of \p circ is applied to 
     *  the qubit 3 of the current circuit, while qubit 1 of \p circ is applied to
     *  the qubit 2 of the current circuit.
     * @param box a ghost gate to indicate boundaries where 
     *  gates from the circuit \p circ are placed into the current circuit.
     * @param flag_inv if true, first get conjugate transpose gates from \p circ. 
     */
    void copy_gates_from(YCCQ circ, YCVI regs_new, 
        YCCB box = std::shared_ptr<const Box__>(nullptr),
        YCB flag_inv = false);

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
     * @param[in] name name of the register;
     * @param[in] n_qubits number of qubits in the register;
     * @param[in] flag_ancilla if true, it is an ancilla register;
     * @return qubits where the register is placed.
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

    /**
     * @brief Set initial amplitudes to specified qubits.
     * @param[in] qb first qubit to set;
     * @param[in] nq number of qubits to set;
     * @param[in] ampl_vec_real vector with real parts of amplitudes of size 2^nq;
     * @param[in] ampl_vec_imag vector with imaginary parts of amplitudes of size 2^nq;
     */
    void set_init_vector(YCI qb, YCI nq, YVQ ampl_vec_real, YVQ ampl_vec_imag);
    void set_init_vector(INIT_STATE__& init_state);
    void reset_init_vector(INIT_STATE__& state);

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
    void read_structure_gate_swap(YISS istr, YCS path_in, YCB flag_inv=false);
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

        oo_layers_->add_gate(oo);

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

        oo_layers_->add_gate(oo);

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

        oo_layers_->add_gate(oo);

        return get_the_circuit();
    }

    /** Set Pauli X gate at \p t target qubit. 
     * CF: X target
     * @param[in] t target qubit;
     * @param[in] cs control qubits;
     * @return pointer to the circuit.
     * */
    inline YQCP x(YCI t, YVIv cs = {}){ return add_sqg<X__>(t, cs); }
    YQCP x(YCVI ts, YVIv cs = {});

    /** Set Pauli Y gate at \p t target qubit. 
     * CF: Y target
     * @param[in] t target qubit;
     * @param[in] cs control qubits;
     * @return pointer to the circuit.
     * */
    inline YQCP y(YCI t, YVIv cs = {}){ return add_sqg<Y__>(t, cs); }
    YQCP y(YCVI ts, YVIv cs = {});

    /** Set Hadamard gate at \p t target qubit.  
     * CF: H target
     * @param[in] t target qubit;
     * @param[in] cs control qubits;
     * @return pointer to the circuit. 
     */
    inline YQCP h(YCI t, YVIv cs = {}){ return add_sqg<H__>(t, cs); }
    YQCP h(YCVI ts, YVIv cs = {});

    /** Set Pauli Z gate at \p t target qubit. 
     * CF: Z target
     * @param[in] t target qubit;
     * @param[in] cs control qubits;
     * @return pointer to the circuit.
     * */
    inline YQCP z(YCI t, YVIv cs = {}){ return add_sqg<Z__>(t, cs); }
    YQCP z(YCVI ts, YVIv cs = {});

    /** Set a Rx-rotation gate.   
     * CF: Rx target angle
     * @param[in] t target qubit;
     * @param[in] a angle to rotate on;
     * @return pointer to the circuit.
     * */
    inline YQCP rx(YCI t, YCQR a, YVIv cs = {}, YCB flag_inv = false){ return add_sq_rg<Rx__>(t, a, cs, flag_inv); }

    /** Set a Ry-rotation gate.   
     * CF: Ry target angle
     * @param[in] t target qubit;
     * @param[in] a angle to rotate on;
     * @return pointer to the circuit.
     * */
    inline YQCP ry(YCI t, YCQR a, YVIv cs = {}, YCB flag_inv = false){ return add_sq_rg<Ry__>(t, a, cs, flag_inv); }

    /** Set a Rz-rotation gate.   
     * CF: Rz target angle
     * @param[in] t target qubit;
     * @param[in] a angle to rotate on;
     * @return pointer to the circuit.
     * */
    inline YQCP rz(YCI t, YCQR a, YVIv cs = {}, YCB flag_inv = false){ return add_sq_rg<Rz__>(t, a, cs, flag_inv); }

    /** Set the Rc-rotation gate: Ry(ay).Rz(az)   
     * CF: Rc target angles
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

    /** @brief Add a swap operator between qubits \p t1 and \p t2. */
    YQCP swap(YCI t1, YCI t2, YVIv cs = {});

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

    /** Circuit wavefunction analysis.
     * @param[in] organize_state input array that indicates how to output every state
     *        e.g. if n_qubits = 4, and organize_state = [2,1,1], then
     *        amplitude  |i3 i2> |i1> |i0>;
     * @param[in] prec input precision of state amplitudes to output;
     * */
    void wavefunction_standard_analysis(
        const std::vector<unsigned>& organize_state, 
        const unsigned& prec = 3
    );

    /**
     * @brief Get strings with all states and states with only nonzero amplitudes;
     * @param[in] organize_state input array that indicates how to output every state
     *        e.g. if n_qubits = 4, and organize_state = [2,1,1], then
     *        amplitude  |i3 i2> |i1> |i0>;
     * @param[in] prec input precision of state amplitudes to output;
     * @param[out] str_wv string with the all states;
     * @param[out] str_wv_nz string only with states that have non-zero amplitudes;
     */
    void get_wavefunction(
        const std::vector<unsigned>& organize_state, 
        YS str_wv, 
        YS str_wv_nz,
        const unsigned& prec = 3
    );

    /**
     * @brief Get a special wavefunction.
     * @param[in] states_to_choose array with bits that a state must have to be chosen:
     *   Elements of this vector must be 0 or 1. 
     *   -1 means that it can be either 0 or 1.
     *   The first element correspond to the most significant qubit.
     * @param[in] organize_state input array that indicates how to output every state
     *        e.g. if n_qubits = 4, and organize_state = [2,1,1], then
     *        amplitude  |i3 i2> |i1> |i0>;
     * @param[out] str_wv_chosen string with chosen states;
     * @param[in] prec input precision of state amplitudes to output;
     */
    void get_special_wavefunction(
        const std::vector<short>& state_to_choose,
        const std::vector<unsigned>& organize_state, 
        YS str_wv_chosen,
        YCU prec = 3 
    );
    void get_special_wavefunction(
        const std::vector<short>& state_to_choose,
        const std::vector<unsigned>& organize_state, 
        YS str_wv_chosen,
        std::list<std::vector<short>>& states_chosen,
        std::vector<Complex>& ampls_chosen,
        YCU prec = 3 
    );

    /**
     * @brief Return the state vector of the circuit. 
     */
    void get_state_vector(YVQ state_real, YVQ state_imag);

    /**
     * @brief Assume that ancillae are the top-priority qubits.
     * Calculates amplitudes of all qubits except the ancillae.
     * @param organize_state input array that indicates how to output every state:
     *        e.g. if n_qubits = 4, and organize_state = [2,1,1], then
     *        amplitude  |i3 i2> |i1> |i0>.
     *        This array must take into account ancillae.
     * @param[out] str_wv_out string with resulting states (set all ancillae to zero);
     * @param[out] states_out resulting states (set all ancillae to zero);
     * @param[out] ampls_out amplitudes of the resulting states;
     * @param[in] state_to_choose array with bits that a state must have to be chosen:
     *   Elements of this vector must be 0 or 1. 
     *   -1 means that it can be either 0 or 1.
     *   The first elements in the array correspond to the high-priority qubits.
     *   The size of the array must take into account the ancillae.
     * @param prec precision of the amplitudes in \p str_wv_out;
     */
    void get_state_zero_ancillae(
        YCVU organize_state, 
        YS str_wv_out,
        std::list<std::vector<short>>& states_out,
        std::vector<Complex>& ampls_out,
        YCVsh state_to_choose = YVshv{},
        YCU prec = 3
    );
    void get_state_full(
        YCVU organize_state, 
        YS str_wv_out,
        std::list<std::vector<short>>& states_out,
        std::vector<Complex>& ampls_out,
        YCVsh state_to_choose = YVshv{},
        YCU prec = 3
    );

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


private:
    // if one addes a new property, do not forget to add it to the copy constructor.

    std::string name_; // name of the circuit;
    QuESTEnv env_; // QuEST execution environment;
    Qureg c_; // quantum register (circuit);
    YMIX::YTimer timer_; // timer for the circuit object;
    YGlobalVariables gv; // global variables;
    std::string path_to_output_; // path where output files to be written
    unsigned nq_; // number of qubits in the circuit;

    std::vector<INIT_STATE__> init_state_; // initial state (as a vector)

    std::string cfname_;  // name of the .circuit file;
    std::string texname_; // name of the .tex file;

    // registers:
    // regs_[rname][i] is the i-th qubit in the register "rname";
    // the 0-th qubit is the least signficant in the register.
    std::map<std::string, YVIv> regs_; 

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
};






#endif