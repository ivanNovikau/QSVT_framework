#include "../include/oracletool.h"
using namespace std;

OracleTool__::OracleTool__(
    const QuESTEnv& env, 
    YCS pname, 
    YCS path_to_inputs
) : BaseTool__(
    env, pname, path_to_inputs, 
    false,  // flag_compute_output (not used here),
    false,  // flag_print_output (not used here),
    false,  // default flag_circuit, 
    false,  // default flag_tex, 
    false,  // default flag_layers, 
    true    // default flag_hdf5
){
    format_file_ = FORMAT_ORACLE;
    sel_compute_output_   = "zero-ancillae";
    sel_print_output_     = "none";
    flag_init_state_file_ = false;

    flag_prob_ = false;

    read_data();
}


OracleTool__::~OracleTool__()
{
    YMIX::print_log( "*** Destruction of the oracle tool is done. ***");
}


void OracleTool__::read_circuit_structure_from_file(YCS data)
{
    istringstream istr(data);
    string word;

    // read data
    while(istr >> word)
    {
        if(YMIX::compare_strings(word, "CONSTANTS"))
            read_constants(istr);

        if(YMIX::compare_strings(word, "OPTIONS"))
            read_options(istr);

        if(YMIX::compare_strings(word, "CIRCUITS_DECLARATION"))
            read_circuit_declaration(istr);

        if(YMIX::compare_strings(word, "CIRCUIT_STRUCTURE"))
            read_circuit_structure(istr);

        if(YMIX::compare_strings(word, "INPUT_STATES"))
            read_input_states(istr);
    }

    YMIX::print_log( "Finish reading file\n");

    // check if the circuit to launch is defined:
    if(!oc_to_launch_)
        throw "Error: The circuit to launch is not defined (there is not a section INPUT_STATES)."s;

    // write circuits to corresponding files:
    for(auto const& it: ocs_)
        it.second->print_gates();

    // set the output format for the circuit of interest:
    oc_to_launch_->set_standart_output_format();

    YMIX::print_log( "Oracle prepared.\n");
}


void OracleTool__::read_constants(YISS istr)
{
    string word, constant_name;
    qreal constant_value;
    YMIX::print_log( "Reading constants...");

    try
    {
        while(istr >> word)
        {
            if(YMIX::compare_strings(word, "END_CONSTANTS"))
                break;

            // name of the constant
            constant_name = word; 

            // value of the constant
            istr >> word;
            constant_value = get_value_from_word(word);

            // save the constant:
            constants_[constant_name] = constant_value;
        }
    }
    catch(YCS e)
    {
        throw "Error in the CONSTANTS section when one reads the constant "s + constant_name + ": "s + e;
    }
}


void OracleTool__::read_options(YISS istr)
{
    string word;
    YMIX::print_log( "Reading options...");

    try
    {
        while(istr >> word)
        {
            if(YMIX::compare_strings(word, "END_OPTIONS"))
                break;

            if(YMIX::compare_strings(word, "sel_compute_output"))
            {
                istr >> sel_compute_output_;
                if(!YMIX::compare_strings(sel_compute_output_, YVSv{"none", "all", "zero-ancillae"}))
                    throw string(
                        "unknown option ["s + sel_compute_output_ + "] for the selector sel_compute_output."s
                    );
                continue;
            }
            if(YMIX::compare_strings(word, "sel_print_output"))
            {
                istr >> sel_print_output_;
                if(!YMIX::compare_strings(sel_print_output_, YVSv{"none", "all", "zero-ancillae"}))
                    throw string(
                        "unknown option ["s + sel_print_output_ + "] for the selector sel_print_output."s
                    );
                continue;
            }

            if(YMIX::compare_strings(word, "flag_circuit"))
            {
                istr >> flag_circuit_;
                continue;
            }

            if(YMIX::compare_strings(word, "flag_tex"))
            {
                istr >> flag_tex_;
                continue;
            }

            if(YMIX::compare_strings(word, "tex_CL"))
            {
                istr >> YGV::tex_circuit_length;
                continue;
            }

            if(YMIX::compare_strings(word, "compute_prob"))
            {
                istr >> flag_prob_;
                if(flag_prob_)
                {
                    oc_to_launch_->read_reg_int(istr, focus_qubits_);
                }
                continue;
            }
        }

        // correct the printing options:
        string sel_comp = "none";
        if(YMIX::compare_strings(sel_compute_output_, sel_comp))
            if(!YMIX::compare_strings(sel_print_output_, sel_comp))
                sel_print_output_ = sel_comp;

        sel_comp = "zero-ancillae";
        if(YMIX::compare_strings(sel_compute_output_, sel_comp))
            if(YMIX::compare_strings(sel_print_output_, "all"))
                sel_print_output_ = sel_comp; // "all" -> "zero-ancillae"

        // output some options:
        YMIX::print_log( "\n--- Initial OPTIONS ---");
        if(!YMIX::compare_strings(sel_compute_output_, "none"))
            YMIX::print_log("-> do not compute output states.");
        if(YMIX::compare_strings(sel_compute_output_, "all"))
            YMIX::print_log("-> compute all output states.");
        if(YMIX::compare_strings(sel_compute_output_, "zero-ancillae are in the zero state."))
            YMIX::print_log("-> compute only output states, where all ancillae are in the zero state.");
        if(flag_prob_)
        {
            stringstream sstr;
            copy(
                focus_qubits_.begin(), 
                focus_qubits_.end(), 
                std::ostream_iterator<int>(sstr, " ")
            );
            YMIX::print_log("-> compute probabilites of the states on the following qubits: " + sstr.str());
        }
            

        if(!flag_circuit_)        YMIX::print_log("-> do not write the " + FORMAT_CIRCUIT + " files.");
        if(!flag_tex_)            YMIX::print_log("-> do not write the " + FORMAT_TEX + " file.");
        if(flag_tex_) 
            YMIX::print_log(
                "-> length of a single row in the .tex file = " + to_string(YGV::tex_circuit_length)
            );
        YMIX::print_log( "---\n");
    }
    catch(YCS e)
    {
        throw "Error in the OPTIONS section: "s + e;
    }
}


qreal OracleTool__::get_value_from_word(YCS word)
{
    if(word.find("<") == string::npos)
    {
        istringstream sstr(word);
        qreal res_value;
        if(!(sstr >> res_value))
            throw "wrong format"s;
        return res_value;
    }

    unsigned first = word.find("<");
    unsigned last = word.find(">");
    string const_name = word.substr(first+1,last-first-1);

    if(constants_.find(const_name) == constants_.end())
        throw "A constant with a name "s + const_name + " is not found."s;

    return constants_[const_name];
}


void OracleTool__::read_circuit_declaration(YISS istr)
{
    string word, circ_name, current_field;
    int n_regs, n_qubits_in_reg;

    YMIX::print_log( "Reading declaration of circuits...");

    try
    {
        while(istr >> word)
        {
            if(YMIX::compare_strings(word, "END_CIRCUITS_DECLARATION"))
                break;

            // name of the circuit
            circ_name = word; 
            
            // number of registers in the circuit
            current_field = "the number of registers"s; 
            istr >> word;
            n_regs = get_value_from_word(word);

            // registers are organized from the high to low priority
            // (from top to bottom qubits):
            vector<int> saved_nqs_in_regs(n_regs);
            vector<string> saved_reg_names(n_regs);
            vector<bool> anc_flags(n_regs);
            int nq_circ = 0;
            YVIv ancs;
            for(unsigned i = 0; i < n_regs; i++)
            {
                // register name
                current_field = "the name of the " + to_string(i) + "-th register";
                istr >> saved_reg_names[i];

                // number of qubits in the register
                current_field = "the number of qubits in the " + to_string(i) + "-th register";
                istr >> word;
                saved_nqs_in_regs[i] = get_value_from_word(word);

                nq_circ += saved_nqs_in_regs[i];

                // ancilla flag
                current_field = "the ancilla flag of the " + to_string(i) + "-th register";
                istr >> word;
                anc_flags[i] = int(get_value_from_word(word));
            }

            // create a circuit if necessary; choose a circuit:
            if(ocs_.find(circ_name) != ocs_.end())
            {
                ostringstream inf;
                inf << "Warning: The circuit " << circ_name <<  " is declared several times. " <<
                        "The first declaration is taken.\n";
                YMIX::print_log( inf.str());          
            }
            else
                ocs_[circ_name] = make_shared<QCircuit>(
                    circ_name, env_, path_inputs_, nq_circ, 
                    constants_,
                    flag_circuit_, flag_tex_, flag_layers_
                );
            YSQ oc = ocs_[circ_name];

            // add the registers to the chosen circuit
            for(unsigned i = 0; i < n_regs; i++)
                oc->add_register(saved_reg_names[i], saved_nqs_in_regs[i], anc_flags[i]);
            oc->save_regs();
        }
    }
    catch(YCS e)
    {
        throw "Error when one reads "s + current_field + " in the circuit "s + circ_name + ": "s + e;
    }

    // information about the declared circuits
    for(const auto& it: ocs_) 
    {
        ostringstream inf;
        inf << it.first << " with " << it.second->get_n_qubits() << " qubit(s):";
        YMIX::print_log( inf.str(), 1); 

        inf.str(""); inf.clear();
        auto reg_names = it.second->get_reg_names();
        inf << "with " << reg_names.size() << " register(s):\n";
        for(auto const& reg_name: reg_names)
            inf << "register " << reg_name << " with " << 
                it.second->get_nq_in_reg(reg_name) << " qubit(s);\n";
        YMIX::print_log( inf.str(), 2); 
    }
}


void OracleTool__::read_circuit_structure(YISS istr)
{
    string word;
    string curr_circuit_name;
    int n_regs, n_qubits_in_reg;

    // choose the circuit according to its name:
    istr >> word;
    YMIX::print_log( "Reading structure of the circuit " + word);
    if(ocs_.find(word) == ocs_.end()) 
    {
        YMIX::print_log(
            "Warning: No circuit with the name " + word + " has been declared. Skip it.", 1
        );
        return;
    }
    curr_circuit_name = word;
    YSQ oc = ocs_[curr_circuit_name];

    // --- gates in the circuit ---
    while(istr >> word)
    {
        if(YMIX::compare_strings(word, "END_CIRCUIT_STRUCTURE"))
            break;

        if(YMIX::compare_strings(word, "gate"))
            read_gate(istr, oc);

        if(YMIX::compare_strings(word, "igate"))
            read_gate(istr, oc, true);

        if(YMIX::compare_strings(word, "circuit"))
            read_subcircuit(istr, oc);

        if(YMIX::compare_strings(word, "icircuit"))
            read_subcircuit(istr, oc, true);
    }
}


void OracleTool__::read_gate(YISS istr, YPQC oc, YCB flag_inv)
{
    string gate_name;
    qreal par_gate, par_gate2;
     
    istr >> gate_name;
    try
    {
        if(oc->read_structure<X__>(gate_name, istr, flag_inv)) return;
        if(oc->read_structure<Y__>(gate_name, istr, flag_inv)) return;
        if(oc->read_structure<Z__>(gate_name, istr, flag_inv)) return;
        if(oc->read_structure<H__>(gate_name, istr, flag_inv)) return;

        if(oc->read_structure<Rx__>(gate_name, istr, par_gate, flag_inv)) return;
        if(oc->read_structure<Ry__>(gate_name, istr, par_gate, flag_inv)) return;
        if(oc->read_structure<Rz__>(gate_name, istr, par_gate, flag_inv)) return;
        if(oc->read_structure<Phase__>(gate_name, istr, par_gate, flag_inv)) return;

        if(oc->read_structure<Rc__>(gate_name, istr, par_gate, par_gate2, flag_inv)) return;


        if(YMIX::compare_strings(gate_name, "condR"))
        {
            oc->read_structure_gate_condR_split(istr, path_inputs_, flag_inv);
            return;
        }
        if(YMIX::compare_strings(gate_name, "adder1"))
        {
            oc->read_structure_gate_adder1(istr, path_inputs_, flag_inv);
            return;
        }
        if(YMIX::compare_strings(gate_name, "adder"))
        {
            oc->read_structure_gate_adder(istr, path_inputs_, flag_inv);
            return;
        }
        if(YMIX::compare_strings(gate_name, "subtractor1"))
        {
            oc->read_structure_gate_subtractor1(istr, path_inputs_, flag_inv);
            return;
        }
        if(YMIX::compare_strings(gate_name, "swap"))
        {
            oc->read_structure_gate_swap(istr, path_inputs_, flag_inv);
            return;
        }
        if(YMIX::compare_strings(gate_name, "Fourier"))
        {
            oc->read_structure_gate_fourier(istr, path_inputs_, flag_inv);
        }
        if(YMIX::compare_strings(gate_name, "PE"))
        {
            oc->read_structure_gate_phase_estimation(istr, path_inputs_, ocs_, flag_inv);
        }

    }
    catch(YCS e)
    {
        ostringstream ostr;
        ostr << "--- Error in the structure of the circuit " << oc->get_name() 
             << ", in the gate " << gate_name << " ---\n" <<
                "" << e << " ---\n";
        throw ostr.str();
    }
}


void OracleTool__::read_subcircuit(YISS istr, YPQC oc, YCB flag_inv)
{
    string subcircuit_name, word;
    string curr_circuit_name = oc->get_name();
    int id_qubit;
    unsigned nq_sub;
    bool flag_skip = false;
    bool flag_plain = false;
    bool flag_end_circuit = false;

    istr >> subcircuit_name;
    if(ocs_.find(subcircuit_name) == ocs_.end())
    {
        string warn_line;
        warn_line = "\n\n-------------------------------------------------------------------------------\n";
        warn_line += "--- Warning: setting the structure of the circuit [" + curr_circuit_name + "] ---\n";
        warn_line += "The subcircuit [" + subcircuit_name + "] is not found. We skip it.\n";
        warn_line += "-------------------------------------------------------------------------------\n";
        YMIX::print_log( warn_line);
        flag_skip = true;
    }
    if(YMIX::compare_strings(subcircuit_name, curr_circuit_name))
    {
        string warn_line;
        warn_line = "\n\n-------------------------------------------------------------------------------\n";
        warn_line += "--- Warning: setting the structure of the circuit [" + curr_circuit_name + "] ---\n";
        warn_line += "Self-insertion: circuit cannot include itself as a subcircuit. We skip it.\n";
        warn_line += "-------------------------------------------------------------------------------\n";
        YMIX::print_log( warn_line);
        flag_skip = true;
    }
    if(ocs_[subcircuit_name]->get_n_gates() == 0)
    {
        string warn_line;
        warn_line = "\n\n-------------------------------------------------------------------------------\n";
        warn_line += "--- Warning: setting the structure of the circuit [" + curr_circuit_name + "] ---\n";
        warn_line += "The subcircuit [" + subcircuit_name + "] is empty. We skip it.\n";
        warn_line += "-------------------------------------------------------------------------------\n";
        YMIX::print_log( warn_line);
        flag_skip = true;
    }

    if(!flag_skip)
    {
        YSQ oc_sub  = ocs_[subcircuit_name];
        int plain_int;

        istr >> word;

        // check whether subcircuit qubits are connected to 
        //         the current circuit qubits in a trivial way:
        try 
        {
            plain_int = stoi(word);
            if(plain_int == -1) flag_plain = true;
        }
        catch(...)
        {   
            flag_plain = false;
        }

        nq_sub = oc_sub->get_n_qubits();
        vector<int> ids_q(nq_sub);
        if(flag_plain)
            for(unsigned i = 0; i < nq_sub; i++)
                ids_q[i] = i;
        else
        {
            ids_q = YVIv {};
            oc->read_reg_int(istr, ids_q, word);
            if(ids_q.size() != nq_sub)
            {
                string err_line;
                err_line  = "--- Error: setting the structure of the circuit " + curr_circuit_name + " ---\n";
                err_line += "The subcircuit " + subcircuit_name + " has " + to_string(nq_sub) + " qubits, while" + 
                    " one has indicated " + to_string(ids_q.size()) + " qubits to connect to.";
                throw err_line;
            }
        }
        oc->copy_gates_from(oc_sub, ids_q, YSB(nullptr), flag_inv);

        // read flag whether it is necessary to show output after the subcircuit:
        int flag_output = 1;
        istr >> word;
        if(YMIX::compare_strings(word, "end_circuit"))
            flag_end_circuit = true;
        else
        {
            try
            {
                flag_output = get_value_from_word(word);
            }  
            catch(...)
            {
                throw "error while reading flag_output for the subcircuit " + subcircuit_name + 
                      ": the word \"" + word + "\" cannot be converted into flag_output.";
            }
        }
            
        // add a Stop gate to output a state after this subcircuit:
        if(flag_output)
        {
            oc->add_stop_gate(oc_sub->get_name());
        } 
    }

    if(!flag_end_circuit)
        while(istr >> word)
            if(YMIX::compare_strings(word, "end_circuit"))
                break;
}


void OracleTool__::read_input_states(YISS istr)
{
    YMIX::print_log( "Reading input states...");
    string word, name_of_main_circuit;

    // define a circuit to launch:
    istr >> name_of_main_circuit;
    if(ocs_.find(name_of_main_circuit) == ocs_.end())
        throw "Error: The circuit with a name \"" + name_of_main_circuit + "\" is not found.\n" + 
            "This circuit cannot be set as a circuit to launch.";
    oc_to_launch_ = ocs_[name_of_main_circuit];
    
    try
    {
        // to read the initial state directly from the .oracle file
        while(istr >> word)
        {
            // to read the initial state from the .init file:
            if(YMIX::compare_strings(word, "INIT_FILE"))
            {
                read_state_init_file();
                break;
            }
                
            if(YMIX::compare_strings(word, "STATE"))
                read_state(istr);

            if(YMIX::compare_strings(word, "END_INPUT_STATES"))
                break;
        }
    }
    catch(YCS e)
    {
        ostringstream ostr;
        ostr << "--- Error while reading initial states: ---"  <<
                "\n--- " << e << " ---\n";
        throw ostr.str();
    }
}


void OracleTool__::read_state(YISS istr)
{
    string word;
    std::map<string, vector<int>> one_state;
    vector<int> ids_qs;

    oc_to_launch_->read_reg_int(istr, ids_qs);
    one_state[YGV::reg_whole_circuit] = ids_qs;
    init_states_.push_back(one_state); 
}


void OracleTool__::read_state_init_file()
{
    int N;
    flag_init_state_file_ = true;

    string fname_init = path_inputs_ + "/" + pname_ + FORMAT_INIT;
    cout << "Reading the file: " << path_inputs_ + "/" + pname_ + FORMAT_INIT << endl;
    YMIX::read_init_state(fname_init, init_ampl_vec_real_, init_ampl_vec_imag_);
}


void OracleTool__::launch()
{
    YMIX::YTimer timer_comp;
    YMIX::print_log( 
        "--- Analysis of the circuit [" + oc_to_launch_->get_name() + "] ---"
    );

    // working circuit object:
    shared_ptr<QCircuit> u_work = oc_to_launch_;

    // ---- store basic data:
    hfo_.open_w();

    // number of qubits
    hfo_.add_scalar(u_work->get_n_qubits(), "nq", "basic");
    hfo_.add_scalar(u_work->get_na(), "na", "basic");

    // register names
    string res_lin = "";
    for(auto const& reg_name: u_work->get_reg_names())
        res_lin += reg_name + ", ";
    res_lin.pop_back(); res_lin.pop_back();
    hfo_.add_scalar(res_lin, "register-names", "basic");

    // number of qubits in every register:
    hfo_.add_vector(u_work->get_standart_output_format(), "register-nq", "basic");

    // number of initial states:
    uint32_t n_init_states = flag_init_state_file_ ? 1: init_states_.size();
    hfo_.add_scalar(n_init_states, "n-init-states", "states");

    // constants:
    for(auto const& [key, val] : constants_)
    {
        hfo_.add_scalar(val, key, "constants");
    }

    hfo_.close();

    // --- Analyse and Store if necessary initial and output states of the circuit ---
    int count_init_state = 0;
    if(flag_init_state_file_)
    {
        // remove gates from the previous launch:
        u_work->reset_qureg();

        cout << "size of the initial vector is " << init_ampl_vec_real_.size() << endl;
        u_work->set_init_vector(init_ampl_vec_real_, init_ampl_vec_imag_);
        calc(u_work, count_init_state, timer_comp);
    }
    else
    {
        for(auto const& state: init_states_)
        {
            // empty previous initial binary states:
            u_work->empty_binary_states();

            // remove gates from the previous launch:
            u_work->reset_qureg();
            
            // set initial states:
            for(auto const& reg: state)
                u_work->set_reg_state(reg.first, reg.second);
            u_work->set_init_binary_state();

            // print input and output states:
            calc(u_work, count_init_state, timer_comp);
            ++count_init_state;
        }
    }
}


void OracleTool__::calc(shared_ptr<QCircuit>& u_work, YCI count_init_state, YMIX::YTimer& timer_comp)
{
    // --- Print the initial state ---
    {
        YMIX::StateVectorOut outF;
        u_work->get_state(outF);
        YMIX::print_log(
            ".....................\n...Initial state " + to_string(count_init_state) + "...\n" + outF.str_wv
        );
    
        // --- Store the initial state ---
        hfo_.open_w();
        hfo_.add_vector(outF.ampls,  "initial-amplitudes-"s + to_string(count_init_state), "states");
        hfo_.add_matrix(outF.states, "initial-states-"s + to_string(count_init_state),     "states");
        hfo_.close(); 
    }

    // --- Print output states ---
    if(!YMIX::compare_strings(sel_compute_output_, "none"))
    {
        int id_current_gate = 0;
        string stop_point_name;
        YMIX::StateVectorOut outF, outZ;
        while(id_current_gate < u_work->get_n_gates())
        {
            // generate the circuit:
            timer_comp.Start();
            YMIX::print_log( "Calculating the circuit... ", 0, false, false);
            u_work->generate(stop_point_name, id_current_gate);

            u_work->get_state(outZ, true);
            if(YMIX::compare_strings(sel_compute_output_, "all"))
                u_work->get_state(outF);

            timer_comp.Stop();
            YMIX::print_log( "duration: " + timer_comp.get_dur_str_s());

            if(YMIX::compare_strings(sel_print_output_, "all"))
                YMIX::print_log(
                        "...Output all states after " + stop_point_name + ": \n" + outF.str_wv
                );
            else if(YMIX::compare_strings(sel_print_output_, "zero-ancillae"))
                YMIX::print_log(
                        "...Output zero-ancilla states after " + stop_point_name + ": \n" + outZ.str_wv
                );
        }

        // --- Store the output state at the very end of the circuit ---
        hfo_.open_w();
        if(YMIX::compare_strings(sel_compute_output_, "all"))
        {
            hfo_.add_vector(
                outF.ampls,  
                "output-all-amplitudes-"s + to_string(count_init_state), 
                "states"
            );
            hfo_.add_matrix(
                outF.states, 
                "output-all-states-"s + to_string(count_init_state),     
                "states"
            );
        }
        if(YMIX::compare_strings(sel_compute_output_, YVSv{"zero-ancillae", "all"}))
            if(outZ.ampls.size() > 0)
            {
                hfo_.add_vector(
                    outZ.ampls,  
                    "output-zero-anc-amplitudes-"s + to_string(count_init_state), 
                    "states"
                );
                hfo_.add_matrix(
                    outZ.states, 
                    "output-zero-anc-states-"s + to_string(count_init_state),     
                    "states"
                );
            }
        hfo_.close(); 
    }

    // --- calculate state probabilities on the indicated qubits ---
    {
        if(YMIX::compare_strings(sel_compute_output_, "none"))
        {
            timer_comp.Start();
            YMIX::print_log( "Calculating the circuit... ", 0, false, false);
            u_work->generate();
            timer_comp.Stop();
            YMIX::print_log( "duration: " + timer_comp.get_dur_str_s());
        }// otherwise, the circuit has been already generated;


        hfo_.open_w();
        if(flag_prob_)
        {
            vector<qreal> outProbs(1<<focus_qubits_.size());
            YMIX::print_log( "Calculating the probabilities... \n", 0, false, false);
            calcProbOfAllOutcomes(
                &outProbs[0], 
                oc_to_launch_->get_qureg(), 
                &focus_qubits_[0],
                focus_qubits_.size()
            );

            hfo_.add_group("probabilities");
            hfo_.add_vector(focus_qubits_,  "qubits"s, "probabilities");
            hfo_.add_vector(outProbs,  "probs", "probabilities");
        }
        hfo_.close(); 
    }
}