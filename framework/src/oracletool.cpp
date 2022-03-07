#include "../include/oracletool.h"
using namespace std;

OracleTool__::OracleTool__(
    const QuESTEnv& env, 
    YCS pname, 
    YCS path_to_inputs, 
    const bool& flag_compute_output,
    const bool& flag_test,
    const bool& flag_compute_iterator_output
){
    flag_test_ = flag_test;
    flag_compute_output_ = flag_compute_output;
    flag_compute_iterator_output_ = flag_compute_iterator_output;

    if(flag_compute_iterator_output_) flag_compute_output_ = false;

    env_ = env;

    pname_ = pname;
    path_inputs_ = path_to_inputs;
    ifname_ = path_inputs_ + "/" + pname_ + FORMAT_ORACLETOOL;

    string data;
    char cdata[YSIZE_CHAR_ARRAY];
    if(env_.rank == 0)
    {
        string current_path = filesystem::current_path();
        YMIX::print_log(env_, "Current path: " + current_path);
        YMIX::print_log(env_, "\n*** Reading an INPUT FILE ***");
        read_input_file(data);
    }

    int size_data;
    if(env_.rank == 0) 
    {
        if(YSIZE_CHAR_ARRAY < data.size()) 
            throw "Error: Size of a char array is too small to save the input file."s;

        strcpy(cdata, data.c_str());
        size_data = data.size();
    }

    if(YMPI)
    {
        MPI_Bcast(&size_data, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
        MPI_Bcast(cdata, size_data+1, MPI_CHAR, 0, MPI_COMM_WORLD);
    }
    
    if(env_.rank > 0) data = string(cdata);

    read_circuit_structure_from_file(data);
}

OracleTool__::~OracleTool__()
{
    YMIX::print_log(env_, "oracletool destruction");
}

void OracleTool__::read_input_file(YS data)
{
    ifstream ff(ifname_);
    if(!ff.is_open()) throw "Error: Here there is not a file: " + ifname_;
    data = string((istreambuf_iterator<char>(ff)), istreambuf_iterator<char>());
    ff.close();

    // clean the buffer from empty lines and comments
    istringstream istr(data);
    string data_clr = "";
    string line;
    while(getline(istr, line))
    {
        line = YMIX::remove_comment(line);
        line = YMIX::trim(line);
        if(line.find_first_not_of(' ') == string::npos)
            continue;
        data_clr += line + "\n";
    }
    std::transform(data_clr.begin(), data_clr.end(), data_clr.begin(), ::tolower);
    data = data_clr;
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

        if(YMIX::compare_strings(word, "CIRCUITS_DECLARATION"))
            read_circuit_declaration(istr);

        if(YMIX::compare_strings(word, "CIRCUIT_STRUCTURE"))
            read_circuit_structure(istr);

        if(YMIX::compare_strings(word, "INPUT_STATES"))
            read_input_states(istr);
    }

    YMIX::print_log(env_, "Finish reading file\n");

    // check if a circuit to launch is defined:
    if(!oc_to_launch_)
        throw "Error: A circuit to launch is not defined (there is not a section INPUT_STATES)."s;

    // write circuits to corresponding files:
    for(auto const& it: ocs_)
        it.second->print_gates(true);

    // set an output format for the circuit of interest:
    oc_to_launch_->set_standart_output_format();

    YMIX::print_log(env_, "Oracle prepared.\n");
}

void OracleTool__::read_constants(YISS istr)
{
    string word, constant_name;
    qreal constant_value;
    YMIX::print_log(env_, "Reading constants...");

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

qreal OracleTool__::get_value_from_word(YCS word)
{
    if(word.find("<") == string::npos)
    {
        istringstream sstr(word);
        qreal res_value;
        sstr >> res_value;
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
    string word, circ_name;
    
    int n_regs, n_qubits_in_reg;

    YMIX::print_log(env_, "Reading declaration of circuits...");
    while(istr >> word)
    {
        if(YMIX::compare_strings(word, "END_CIRCUITS_DECLARATION"))
            break;

        // name of the circuit
        circ_name = word; 
        
        // number of registers in the circuit
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
            istr >> saved_reg_names[i];

            // number of qubits in the register
            istr >> word;
            saved_nqs_in_regs[i] = get_value_from_word(word);

            nq_circ += saved_nqs_in_regs[i];

            // ancilla flag
            istr >> word;
            anc_flags[i] = int(get_value_from_word(word));
        }

        // create a circuit if necessary; choose a circuit:
        if(ocs_.find(circ_name) != ocs_.end())
        {
            ostringstream inf;
            inf << "Warning: The circuit " << circ_name <<  " is declared several times. " <<
                    "The first declaration is taken.\n";
            YMIX::print_log(env_, inf.str());          
        }
        else
            ocs_[circ_name] = make_shared<QCircuit>(circ_name, env_, path_inputs_, nq_circ, constants_);
        YSQ oc = ocs_[circ_name];

        // add the registers to the chosen circuit
        for(unsigned i = 0; i < n_regs; i++)
            oc->add_register(saved_reg_names[i], saved_nqs_in_regs[i], anc_flags[i]);
        oc->save_reg_names();
    }

    // information about the declared circuits
    for(const auto& it: ocs_) 
    {
        ostringstream inf;
        inf << it.first << " with " << it.second->get_n_qubits() << " qubit(s):";
        YMIX::print_log(env_, inf.str(), 1); 

        inf.str(""); inf.clear();
        auto reg_names = it.second->get_reg_names();
        inf << "with " << reg_names.size() << " register(s):\n";
        for(auto const& reg_name: reg_names)
            inf << "register " << reg_name << " with " << 
                it.second->get_nq_in_reg(reg_name) << " qubit(s);\n";
        YMIX::print_log(env_, inf.str(), 2); 
    }
}

void OracleTool__::read_circuit_structure(YISS istr)
{
    string word;
    string curr_circuit_name;
    int n_regs, n_qubits_in_reg;

    // choose the circuit according to its name:
    istr >> word;
    YMIX::print_log(env_, "Reading structure of a circuit " + word);
    if(ocs_.find(word) == ocs_.end()) 
    {
        YMIX::print_log(
            env_, 
            "Warning: No circuit with a name " + word + " has been declared. Skip it.", 1);
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
        if(oc->read_structure<Z__>(gate_name, istr, flag_inv)) return;
        if(oc->read_structure<H__>(gate_name, istr, flag_inv)) return;

        if(oc->read_structure<Ry__>(gate_name, istr, par_gate, flag_inv)) return;
        if(oc->read_structure<Rz__>(gate_name, istr, par_gate, flag_inv)) return;
        if(oc->read_structure<Phase__>(gate_name, istr, par_gate, flag_inv)) return;

        if(oc->read_structure<Rc__>(gate_name, istr, par_gate, par_gate2, flag_inv)) return;


        if(YMIX::compare_strings(gate_name, "condR"))
        {
            //oc->read_structure_gate_condR(istr, path_inputs_, flag_inv);
            oc->read_structure_gate_condR_split(istr, path_inputs_, flag_inv);
            return;
        }
        if(YMIX::compare_strings(gate_name, "adder1"))
        {
            oc->read_structure_gate_adder1(istr, path_inputs_, flag_inv);
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

    }
    catch(YCS e)
    {
        ostringstream ostr;
        ostr << "--- Error in a structure of a circuit " << oc->get_name() 
             << ", in a gate " << gate_name << " ---\n" <<
                "--- " << e << " ---\n";
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
        warn_line  = "--- Warning: setting the structure of a circuit " + curr_circuit_name + " ---\n";
        warn_line += "The subcircuit " + subcircuit_name + " is not found. We skip it.\n";
        warn_line += "--------------------------------------------------\n";
        YMIX::print_log(env_, warn_line);
        flag_skip = true;
    }
    if(YMIX::compare_strings(subcircuit_name, curr_circuit_name))
    {
        string warn_line;
        warn_line  = "--- Warning: setting the structure of a circuit " + curr_circuit_name + " ---\n";
        warn_line += "Self-insertion: the circuit has itself as a subcircuit. We skip it.\n";
        warn_line += "--------------------------------------------------\n";
        YMIX::print_log(env_, warn_line);
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
                err_line  = "--- Error: setting the structure of a circuit " + curr_circuit_name + " ---\n";
                err_line += "The subcircuit " + subcircuit_name + " has " + to_string(nq_sub) + " qubits, while" + 
                    " one has indicated " + to_string(ids_q.size()) + " qubits to connect to.";
                throw err_line;
            }
        }
        oc->copy_gates_from(oc_sub, ids_q, YSB(nullptr), flag_inv);

        // read flag whether it is necesasry to show output after the subcircuit:
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
                throw "error while reading flag_output for a subcircuit " + subcircuit_name + 
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
    YMIX::print_log(env_, "Reading input states...");
    string word, name_of_main_circuit;

    // define a circuit to launch:
    istr >> name_of_main_circuit;
    if(ocs_.find(name_of_main_circuit) == ocs_.end())
        throw "Error: The circuit with a name \"" + name_of_main_circuit + "\" is not found.\n" + 
            "This circuit cannot be set as a circuit to launch.";
    oc_to_launch_ = ocs_[name_of_main_circuit];
    
    try
    {
        while(istr >> word)
        {
            if(YMIX::compare_strings(word, "STATE"))
                read_state(istr);

            if(YMIX::compare_strings(word, "END_INPUT_STATES"))
                break;
        }
    }
    catch(const string& e)
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

    if(flag_test_) YMIX::print_log(env_, "State number " + to_string(init_states_.size()));
    oc_to_launch_->read_reg_int(istr, ids_qs);
    one_state[gv.reg_whole_circuit] = ids_qs;
    init_states_.push_back(one_state); 
}

// void OracleTool__::launch()
// {
//     string str_wv, str_wv_nz;
//     YMIX::YTimer timer_gen;
//     YMIX::YTimer timer_comp;

//     YMIX::print_log(
//         env_, 
//         "--- Analysis of a circuit " + oc_to_launch_->get_name() + " ---"
//     );

//     // work circuit object:
//     shared_ptr<QCircuit> u_work;

//     // oracle:
//     if(flag_compute_output_)
//     {
//         u_work = 
//     }

//     // create an iterator:
//     if(flag_compute_iterator_output_)
//     {
//         auto cu = make_shared<QCircuit>(oc_to_launch_);
//         cu->controlled(oc_to_launch_->get_n_qubits());

//     }


//     int count_init_state = 0;
//     for(auto const& state: init_states_)
//     {
//         // empty previous initial binary states:
//         oc_to_launch_->empty_binary_states();

//         // remove gates from the previous launch:
//         oc_to_launch_->reset_qureg();
        
//         // set initial states:
//         for(auto const& reg: state)
//             oc_to_launch_->set_reg_state(reg.first, reg.second);
//         oc_to_launch_->set_init_binary_state();

//         // get 
//         oc_to_launch_->get_wavefunction(
//             oc_to_launch_->get_standart_output_format(), str_wv, str_wv_nz, 3
//         );
//         YMIX::print_log(
//             env_,
//             ".....................\n...Initial state " + to_string(count_init_state) + "...\n" + str_wv_nz
//         );

//         // output from an oracle
//         if(flag_compute_output_)
//         {
//             int id_current_gate = 0;
//             string stop_point_name;
//             while(id_current_gate < oc_to_launch_->get_n_gates())
//             {
//                 // generate the circuit:
//                 timer_gen.Start();
//                 YMIX::print_log(env_, "Circuit generation... ", 0, false, false);
//                 oc_to_launch_->generate(stop_point_name, id_current_gate, false);
//                 timer_gen.Stop();
//                 YMIX::print_log(env_, "duration: " + timer_gen.get_dur_str_s());

//                 // compute the output state
//                 timer_comp.Start();
//                 YMIX::print_log(env_, "Circuit computation... ", 0, false, false);
//                 oc_to_launch_->get_wavefunction(
//                     oc_to_launch_->get_standart_output_format(), str_wv, str_wv_nz, 3
//                 );
//                 timer_comp.Stop();
//                 YMIX::print_log(env_, "duration: " + timer_comp.get_dur_str_s());

//                 YMIX::print_log(env_, "...Output state after " + stop_point_name + ": \n" + str_wv_nz);
//             }
//         }
//         ++count_init_state;
//     }
// }


void OracleTool__::launch()
{
    string str_wv, str_wv_nz;
    YMIX::YTimer timer_gen;
    YMIX::YTimer timer_comp;

    YMIX::print_log(
        env_, 
        "--- Analysis of a circuit " + oc_to_launch_->get_name() + " ---"
    );

    // work circuit object:
    shared_ptr<QCircuit> u_work;
    string str_work;

    // oracle:
    if(flag_compute_output_)
    {
        u_work = oc_to_launch_;
        str_work = "Oracle";
    }

    // create the iterate W:
    if(flag_compute_iterator_output_)
    {
        // qubits from the oracle:
        map<string, unsigned> regs_oracle;
        unsigned nqw = 0;
        for(auto const& reg_name: oc_to_launch_->get_reg_names())
        {
            auto reg_nq = oc_to_launch_->get_nq_in_reg(reg_name);
            regs_oracle[reg_name] = reg_nq;
            nqw += reg_nq;
        }

        // controlled oracle
        auto cu = make_shared<QCircuit>(oc_to_launch_);
        cu->controlled(nqw);

        // iterate object:
        nqw += 1; // include q-register
        u_work = make_shared<QCircuit>("W", env_, path_inputs_, nqw);

        map<string, vector<int>> regs;
        regs["q"] = u_work->add_register("q", 1);
        for(auto const& reg_name: oc_to_launch_->get_reg_names())
            regs[reg_name] = u_work->add_register(reg_name, regs_oracle[reg_name]); 
        u_work->save_reg_names();
        u_work->set_standart_output_format();

        auto ts_box = YMATH::get_range(0,nqw-1);
        YVIv cs_box = regs["q"];
        auto rq = regs["q"][0];
        auto ancs = oc_to_launch_->get_ancillae();

        // gates of the iterate:
        u_work->h(rq); // for outputing only

        // add cU 
        u_work->x(rq);
        u_work->copy_gates_from(cu, YMATH::get_range(0, nqw), make_shared<Box__>("U", ts_box, cs_box));
        u_work->x(rq);

        // add cU*
        auto conjU = make_shared<QCircuit>(cu);
        conjU->conjugate_transpose();
        u_work->copy_gates_from(conjU, YMATH::get_range(0, nqw), make_shared<Box__>("U*", ts_box, cs_box));

        // add a reflector operator 
        // cw_->x(rq);  // operator S !!! QSP works without this operator !!!
        u_work->h(rq);
        u_work->x(rq);

        for(auto& anc1: ancs)
            u_work->x(anc1);
        u_work->phase(rq, M_PI, ancs);
        for(auto& anc1: ancs)
            u_work->x(anc1);
        u_work->phase(rq, M_PI);

        u_work->x(rq);
        u_work->h(rq);


        u_work->h(rq); // for outputing only
        str_work = "Iterate";
    }


    int count_init_state = 0;
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

        // get 
        u_work->get_wavefunction(
            u_work->get_standart_output_format(), str_wv, str_wv_nz, 3
        );
        YMIX::print_log(
            env_,
            ".....................\n...Initial state " + to_string(count_init_state) + "...\n" + str_wv_nz
        );

        // output from an oracle
        if(flag_compute_output_ || flag_compute_iterator_output_)
        {
            int id_current_gate = 0;
            string stop_point_name;
            while(id_current_gate < u_work->get_n_gates())
            {
                // generate the circuit:
                timer_gen.Start();
                YMIX::print_log(env_, str_work + " generation... ", 0, false, false);
                u_work->generate(stop_point_name, id_current_gate, false);
                timer_gen.Stop();
                YMIX::print_log(env_, "duration: " + timer_gen.get_dur_str_s());

                // compute the output state
                timer_comp.Start();
                YMIX::print_log(env_, str_work + " computation... ", 0, false, false);
                u_work->get_wavefunction(
                    u_work->get_standart_output_format(), str_wv, str_wv_nz, 3
                );
                timer_comp.Stop();
                YMIX::print_log(env_, "duration: " + timer_comp.get_dur_str_s());

                YMIX::print_log(env_, "...Output state after " + stop_point_name + ": \n" + str_wv_nz);
            }
        }
        ++count_init_state;
    }
}