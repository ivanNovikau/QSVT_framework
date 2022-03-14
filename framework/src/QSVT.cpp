#include "../include/QSVT.h"

using namespace std;


QSVT__::QSVT__(const QuESTEnv& env, YCS project_name, YCS path_inputs)
 :env_(env), 
 project_name_(project_name),
 fname_input_(path_inputs + "/" + project_name + FORMAT_QSP), 
 fname_init_(path_inputs + "/" + project_name + FORMAT_INIT),
 path_inputs_(filesystem::canonical(path_inputs)),
 flag_restart_(false)
{
    timer_.Start();

    flag_circuit_ = false;
    flag_print_zero_states_ = true;
    flag_print_all_states_ = false;
    flag_layers_ = false;

    // --- Set the name of the restart file ---
    rf_.set_name(path_inputs_ + "/" + project_name_ + ENDING_FORMAT_RESTART);
}


QSVT__::~QSVT__()
{
    YMIX::print_log(env_, "Destruction of the QSVT framework...");

    // destroy inner circuits:
    u_.reset();
    iu_.reset();
    codd_.reset();
    ceven_.reset();
    codd_controlled_.reset();
    ceven_controlled_.reset();

    // calculation time
    timer_.Stop();
    ostringstream ofs;
    ofs << "Total simulation time of the QSVT circuit: " << scientific << setprecision(3) 
        << timer_.get_dur() << " ms; "
        << timer_.get_dur() / (60.e3) << " minutes.";
    YMIX::print_log(env_, ofs.str());

    // save time of calculation to the .hdf5 file:
    stringstream sstr;
    sstr << scientific << setprecision(3) 
        << timer_.get_dur() << " ms; "
        << timer_.get_dur() / (60.e3) << " minutes.";
    hfo_.open_w();
    hfo_.add_scalar(sstr.str(), "total-simulation-time", "basic");
    hfo_.close();
}


void QSVT__::read_main_parameters()
{
    ifstream ff(fname_input_);
    if(!ff.is_open())
        YMIX::print_log_err(env_, "Error: there is not a file: " + fname_input_);

    string line;
    string key_name;
    while (getline(ff, line))
    {
        key_name = "";
        line = YMIX::remove_comment(line);
        if(line.find_first_not_of(' ') == string::npos)
            continue;

        istringstream iss(line);
        iss >> key_name;

        // QSVT circuit precision:
        if(YMIX::compare_strings(key_name, "eps"))
        {
            iss >> eps_;
            continue;
        }

        // selector how to create the initial state:
        if(YMIX::compare_strings(key_name, "sel_init"))
        {
            string temp;
            iss >> temp;
            if(YMIX::compare_strings(temp, "vector"))
                sel_init_ = SEL_INIT_STATE_PREP::use_init_vector;
            if(YMIX::compare_strings(temp, "oracle"))
                sel_init_ = SEL_INIT_STATE_PREP::use_init_oracle;

            iss >> coef_time_norm_;  
            continue;
        }

        // print or not the QSP (QSVT) .circuit file:
        if(YMIX::compare_strings(key_name, "flag_circuit"))
        {
            string temp;
            iss >> temp;
            if(YMIX::compare_strings(temp, "true"))
                flag_circuit_ = true;
            else
                flag_circuit_ = false;
            continue;
        }

        // whether to print zero-ancilla states to console:
        if(YMIX::compare_strings(key_name, "flag_print_zero_states"))
        {
            string temp;
            iss >> temp;
            if(YMIX::compare_strings(temp, "true"))
                flag_print_zero_states_ = true;
            else
                flag_print_zero_states_ = false;
            continue;
        }

        // whether to print all states to console:
        if(YMIX::compare_strings(key_name, "flag_print_all_states"))
        {
            string temp;
            iss >> temp;
            if(YMIX::compare_strings(temp, "true"))
                flag_print_all_states_ = true;
            else
                flag_print_all_states_ = false;
            continue;
        }

        // flag to read restart data:
        if(YMIX::compare_strings(key_name, "flag_restart"))
        {
            iss >> flag_restart_;
            continue;
        }

        // flag to calculate layers:
        if(YMIX::compare_strings(key_name, "flag_layers"))
        {
            string temp;
            iss >> temp;
            if(YMIX::compare_strings(temp, "true"))
                flag_layers_ = true;
            else
                flag_layers_ = false;
            continue;
        }

        if(read_special_parameters(iss, key_name))
            continue;
    }
    ff.close();
}


void QSVT__::prepare_hdf5_files()
{
    // --- Read a restart counter from a restart file ---
    string str_res_counter = "";
    restart_counter_ = 0;
    if(flag_restart_)
    {
        YMIX::print_log(env_, "Initial state is to be read from the restart file.");
        rf_.open_r();
        rf_.read_scalar(restart_counter_, "restart-counter", "parameters");
        str_res_counter = "_" + to_string(restart_counter_);
        rf_.close();
    }
    string fname_output = path_inputs_ + "/" + project_name_ + str_res_counter + ENDING_FORMAT_OUTPUT;

    // --- Create an output .hdf5 file ---
    YMIX::print_log(env_, "Creating the output .hdf5 file...");
    hfo_.create(fname_output);
    
    // add directories to the file:
    hfo_.add_group("basic"); // paths, basic circuit description, qsp parameters
                            // qsp parameters, angles, number of gates
    hfo_.add_group("states"); // initial and output states

    // date of simulation:
    string str_date_time;
    YMIX::get_current_date_time(str_date_time);
    hfo_.add_scalar(str_date_time, "date-of-simulation", "basic");

    // save the project name and the path to input files:
    hfo_.add_scalar(project_name_, "project-name", "basic");
    hfo_.add_scalar(path_inputs_, "path-inputs", "basic");
    hfo_.add_scalar(filesystem::current_path(), "launch-path", "basic");

    // close the file:
    hfo_.close();
}


void QSVT__::print_init_data()
{
    YMIX::LogFile cf;
    cf << "\n--- Global circuit: " << oc_->get_name() << "--- \n";
    cf << "Number of qubits: " << nq_ << "\n";
    cf << "Number of ancillae: " << oc_->get_na() << "\n";
    cf << "Precision: " << eps_ << ";\n";
    cf << "Register positions:\n";
    oc_->print_reg_positions(cf.of);
}


void QSVT__::save_basic_data()
{
    hfo_.open_w();

    // number of qubits
    hfo_.add_scalar(nq_, "nq", "basic");
    hfo_.add_scalar(oc_->get_na(), "na", "basic");

    // register names
    string res_lin = "";
    for(auto const& reg_name: oc_->get_reg_names())
        res_lin += reg_name + ", ";
    res_lin.pop_back(); res_lin.pop_back();
    hfo_.add_scalar(res_lin, "register-names", "basic");

    // number of qubits in every register:
    hfo_.add_vector(oc_->get_standart_output_format(), "register-nq", "basic");

    // error of the polynomial approximation:
    hfo_.add_scalar(eps_, "qsvt-error", "basic");

    save_basic_specific_data();
        
    hfo_.close();
}


void QSVT__::append(QCircuit* c)
{
    oc_->insert_gates_from(c);
}


void QSVT__::read_angles_def_parity(YCS line_parity, YVQ phis, uint32_t& N_angles)
{
    string fname_full, fname_res;
    ifstream ff;

    qreal f_par, eps;

    fname_res = 
        "angles_" + line_parity + "_param" + to_string(int(1000*f_par_)) + 
        "_eps" + to_string(int(-log10(eps_))) + 
        FORMAT_ANGLES;
    fname_full = path_inputs_ + "/" + fname_res;

    YMIX::print_log(env_, "Angles are read from the file:\n"s + fname_full);
    ff.open(fname_full.c_str());
    if(ff.is_open())
    {
        ff >> f_par >> eps >> N_angles;
        phis = YVQv(N_angles);
        for(unsigned i = 0; i < N_angles; ++i)
            ff >> phis[i];
        ff.close();
    }
    else
    {
        YMIX::print_log_err(env_, "There is not the following file\n" + fname_full);
    }

    // --- recheck the parameters ---
    if(!YMATH::is_zero(abs(eps - eps_)))
    {
        throw string("eps does not correspond to the QSVD angles.");
    }
    if(!YMATH::is_zero(abs(f_par - f_par_)))
    {
        throw string("QSVD parameter (time step or kappa or etc.) does not correspond to the QSVD angles.");
    }

    // --- Print parameters ---
    ostringstream ocs;
    ocs << "--- Parity: " << line_parity << " ---\n";
    ocs << "function parameter = " << f_par << ";   ";
    ocs << "eps = "            << eps << ";   ";
    ocs << "N of angles = "    << N_angles << ";\n";
    YMIX::print_log(env_, "\n"s + ocs.str());

    // --- Save the angles to .hdf5 file ---
    YMIX::print_log(env_, "\nWriting angles to .hdf5 file");
    hfo_.open_w();
    hfo_.add_vector(phis, line_parity + "-angles", "basic");
    hfo_.close();
}


void QSVT__::read_block_encoding_oracle()
{
    YMIX::print_log(env_, "Reading the block-encoding oracle...");
    OracleTool__* oo = new OracleTool__(
        env_, 
        project_name_, 
        path_inputs_, 
        false,
        false,
        false,
        false
    );
    u_ = oo->get_oracle_copy();
    delete oo;
}


void QSVT__::create_circuit()
{ 
    auto na_qsvt = nq_; // number of ancillae specific to the QSVT circuit;

    // qubits from the oracle:
    map<string, unsigned> regs_oracle;
    int reg_nq = 0;
    for(auto const& reg_name: u_->get_reg_names())
    {
        reg_nq = u_->get_nq_in_reg(reg_name);
        regs_oracle[reg_name] = reg_nq;
        nq_ += reg_nq;
    }

    // --- create an empty circuit ---
    YMIX::print_log(env_, "Initialize the framework circuit...");
    oc_ = make_unique<QCircuit>(
        type_, env_, path_inputs_, nq_,
        map<string, qreal>(),
        flag_circuit_, false, flag_layers_
    );

    // add registers:
    regs_["qb"] = oc_->add_register("qb", na_qsvt);
    for(auto const& reg_name: u_->get_reg_names())
        regs_[reg_name] = oc_->add_register(reg_name, regs_oracle[reg_name]);
    oc_->save_regs();
    oc_->set_standart_output_format();

    // indicate the position of ancillae:
    auto ancs = u_->get_ancillae();
    ancs.insert(
        ancs.end(), 
        regs_["qb"].begin(), 
        regs_["qb"].end()
    );
    oc_->set_ancillae(ancs);

    // ---Create the complex_conjugated oracle ---
    auto conjU = make_shared<QCircuit>(u_);
    conjU->conjugate_transpose();
    iu_ = make_shared<QCircuit>(conjU);

    // --- Create parts of the QSVT circuit ---
    create_circuit_components();
}


void QSVT__::create_circuit_component_def_parity(
    YCS line_parity, YCVQ phis, YCU N_angles, 
    YSQ& circ,  
    YCB flag_imaginary
){
    YMIX::YTimer timer;
    qreal aa;

    auto a_ancs = u_->get_ancillae();
    auto na_qsvt = regs_["qb"].size();

    string str1 = "-> Circuit component: " + line_parity;
    if(flag_imaginary) str1 += ", imaginary.";
    YMIX::print_log(env_, str1);

    circ = make_shared<QCircuit>(
        "C-" + line_parity, env_, path_inputs_, nq_,
        map<string, qreal>(),
        flag_circuit_, false, flag_layers_
    );
    auto  q      = circ->add_register("a-qsvt", na_qsvt)[0]; // get the least singificant QSVT ancilla;
    auto  ureg   = circ->add_register("ureg", nq_ - na_qsvt); 
    circ->save_regs();

    // --- Structure the circuit component ---
    timer.StartPrint(env_, "creating... ");

    circ->h(q);
    if(flag_imaginary)
        circ->z(q);

    aa = 2*phis[0];
    circ->x(a_ancs);
    circ->x(q, a_ancs)->rz(q, aa)->x(q, a_ancs);
    circ->x(a_ancs);
    for(uint32_t count_angle = 1; count_angle < N_angles; ++count_angle)
    {
        circ->insert_gates_from(
            u_.get(), 
            make_shared<Box__>("U", ureg, YVIv {})
        );
        circ->z(q);

        aa = 2*phis[count_angle];
        circ->x(a_ancs);
        circ->x(q, a_ancs)->rz(q, aa)->x(q, a_ancs);
        circ->x(a_ancs);
        
        count_angle += 1;
        if(count_angle < N_angles)
        {
            circ->insert_gates_from(
                iu_.get(), 
                make_shared<Box__>("iU", ureg, YVIv {})
            );
            circ->z(q);

            aa = 2*phis[count_angle];
            circ->x(a_ancs);
            circ->x(q, a_ancs)->rz(q, aa)->x(q, a_ancs);
            circ->x(a_ancs);
        }
    }
    circ->h(q);
    timer.StopPrint(env_);

    circ->print_gates();
}


void QSVT__::create_controlled_component(const YSQ circ, YSQ& circ_controlled)
{
    auto b = nq_ - 1;
    circ_controlled = make_shared<QCircuit>(
        circ, circ->get_name() + "-CONTROLLED"
    );
    circ_controlled->controlled({b});
    circ_controlled->print_gates();
}


void QSVT__::set_init_vector()
{
    // -------------------------------------------------------------------------
    // --- Initialization via the oracle ---
    // -------------------------------------------------------------------------
    if(!flag_restart_ && sel_init_ == SEL_INIT_STATE_PREP::use_init_oracle)
    {
        YMIX::print_log(env_, "----------------------------------------------------");
        YMIX::print_log(env_, "Reading oracle for the circuit initialization...");
        YMIX::print_log(env_, "----------------------------------------------------");
        OracleTool__* oo = new OracleTool__(
            env_, 
            project_name_ + "_init", 
            path_inputs_, 
            false,
            false,
            false,
            false
        );
        oc_init_ = oo->get_oracle_copy();
        delete oo;

        // the initialization oracle is inserted into the QSVT circuit qubit by qubit 
        //  starting from the bottom (starting from the least significant qubit):
        auto ts_box = YMATH::get_range(0, oc_init_->get_n_qubits());
        oc_->insert_gates_from(
            oc_init_.get(), 
            make_shared<Box__>("INIT", ts_box, YVIv {})
        );
        oc_->generate();
    }

    // -------------------------------------------------------------------------
    // --- Initialization via the initial vector ---
    // -------------------------------------------------------------------------
    if(!flag_restart_ && sel_init_ == SEL_INIT_STATE_PREP::use_init_vector)
    {
        YMIX::print_log(env_, "----------------------------------------------------");
        YMIX::print_log(env_, "Circuit initialization via a pre-calculated vector...");
        YMIX::print_log(env_, "----------------------------------------------------");

        int N, nq;
        vector<qreal> ampl_vec_real;
        vector<qreal> ampl_vec_imag;

        // read text from the file
        string data;
        {
            ifstream ff(fname_init_);
            if(!ff.is_open()) throw "Error: there is not a file: " + fname_init_;
            data = string((istreambuf_iterator<char>(ff)), istreambuf_iterator<char>());
            ff.close();
        }
        
        // clean the buffer from empty lines and comments
        {
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
        
        // read the data
        string key_name;
        istringstream iss(data);
        while(iss >> key_name)
        {
            // time normalization factor
            if(YMIX::compare_strings(key_name, "beta"))
            {
                iss >> coef_time_norm_;
                continue;
            }

            // number of elements in a vector:
            if(YMIX::compare_strings(key_name, "N"))
            {
                iss >> N;
                ampl_vec_real = vector<qreal>(N);
                ampl_vec_imag = vector<qreal>(N);
                continue;
            }

            // real part of an initial state
            if(YMIX::compare_strings(key_name, "real"))
            {
                for(unsigned i = 0; i < N; ++i)
                    iss >> ampl_vec_real[i];
                continue;
            }

            // imaginary part of an initial state
            if(YMIX::compare_strings(key_name, "imag"))
            {
                for(unsigned i = 0; i < N; ++i)
                    iss >> ampl_vec_imag[i];
                continue;
            }
        }

        nq = log2(N);
        oc_->set_init_vector(0, nq, ampl_vec_real, ampl_vec_imag);
    }

    // -------------------------------------------------------------------------
    // --- Initialization via the restart fiel ---
    // -------------------------------------------------------------------------
    if(flag_restart_)
    {
        YMIX::print_log(env_, "----------------------------------------------------");
        YMIX::print_log(env_, "Initial state is read from the restart file...");
        YMIX::print_log(env_, "----------------------------------------------------");

        int N, nq;
        vector<qreal> ampl_vec_real;
        vector<qreal> ampl_vec_imag;

        // --- read the initial state from the restart file ---
        rf_.open_r();
        rf_.read_scalar(coef_time_norm_, "coef-time-norm", "parameters");
        rf_.read_vector(ampl_vec_real,   "real", "state");
        rf_.read_vector(ampl_vec_imag,   "imag", "state");
        rf_.close();

        N = ampl_vec_real.size();
        nq = log2(N);
        oc_->set_init_vector(0, nq, ampl_vec_real, ampl_vec_imag);
    }

    // -------------------------------------------------------------------------
    // --- Save the initial state ---
    // -------------------------------------------------------------------------
    std::string str_out;
    list<vector<short>> states;
    vector<Complex> ampls;
    oc_->get_state_full(
        oc_->get_standart_output_format(), 
        str_out, 
        states,
        ampls,
        YVshv {},
        3
    );
    YMIX::print_log(env_, "Initial state: \n"s + str_out);
    hfo_.open_w();
    hfo_.add_scalar(coef_time_norm_, "normalization-coef", "basic");
    hfo_.add_vector(ampls, "initial-amplitudes", "states");
    hfo_.add_matrix(states, "initial-states", "states");
    hfo_.close(); 
}


void QSVT__::launch()
{
    string line = "\n\n-------------------\n";
    line           += "--- Computation ---\n";
    line           += "-------------------\n";
    YMIX::print_log(env_, line);

    simulation();
    save_restart_data();
}


void QSVT__::save_restart_data()
{
    YMIX::print_log(env_, "Rewrite the restart .hdf5 file...");
    rf_.create(""); 

    rf_.add_group("parameters");
    rf_.add_group("state");
    
    // save the restart counter:
    ++restart_counter_;
    rf_.add_scalar(restart_counter_, "restart-counter", "parameters");

    // save parameters
    rf_.add_scalar(coef_time_norm_, "coef-time-norm", "parameters");

    // save the last state of the circuit
    vector<qreal> state_real, state_imag;
    oc_->get_state_vector(state_real, state_imag);
    rf_.add_vector(state_real, "real", "state");
    rf_.add_vector(state_imag, "imag", "state");

    rf_.close();
}


void QSVT__::print_gates()
{
    oc_->print_gates();
}








