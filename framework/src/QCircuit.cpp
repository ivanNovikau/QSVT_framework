#include "../include/QCircuit.h"

using namespace std;

QCircuit::QCircuit(YCS name, const QuESTEnv& env, YCS path_to_output, 
    YCU nq, const std::map<std::string, qreal>& constants
){
    timer_.Start();
    name_ = name;
    path_to_output_ = path_to_output;
    env_ = env;
    nq_ = 0;
    if(nq > 0)
        create(nq); 
    else   
        nq_ = nq; 

    constants_ = map<string, qreal>(constants);
}


void QCircuit::create(YCU nq)
{
    if(nq <= 0)
    {
        cerr << "Error: circuit " << name_ << ": number of qubits has to be greater than zero" << endl;
        exit(-1);
    }

    if(nq_ > 0)
    {
        cerr << "Error: circuit " << name_ << ": the circuit has been already created" << endl;
        exit(-1);
    }

    nq_ = nq;
    c_ = createQureg(nq_, env_);

    tex_noc_.resize(nq_);
    tex_lines_.resize(nq_);

    create_circ_file();
    create_tex_file();
    initZeroState(c_);

    oo_layers_ = make_shared<CircuitLayers__>(nq_);
}


QCircuit::QCircuit(YCCQ oc, YCS cname)
{
    timer_.Start();

    if(cname.empty())
        name_ = oc->name_ + "-copy";
    else
        name_ = cname;
    env_ = oc->env_;
    nq_ = oc->nq_;
    c_ = createQureg(nq_, env_);
    initZeroState(c_);

    path_to_output_ = oc->path_to_output_;
    init_state_     = vector<INIT_STATE__>(oc->init_state_);
    regs_           = map<string, YVIv>(oc->regs_);
    regnames_       = vector<string>(oc->regnames_);
    ancs_           = YVIv(oc->ancs_);
    ib_state_       = vector<short>(oc->ib_state_);
    id_start_       = oc->id_start_;
    standart_output_format_ = oc->standart_output_format_;

    tex_noc_ = vector<uint64_t>(tex_noc_);
    tex_lines_ = vector<vector<string>>(tex_lines_);

    create_circ_file();
    create_tex_file();
    save_regs();

    for(auto const& gate: oc->gates_)
    {
        auto gate_copy = gate->copy_gate();
        gates_.push_back(gate_copy);
    }
    constants_ = map<string, qreal>(oc->constants_);
    oo_layers_ = make_shared<CircuitLayers__>(oc->oo_layers_);
}


QCircuit::~QCircuit()
{
    finish_tex_file();
    destroyQureg(c_, env_);
    timer_.Stop();
    if(env_.rank == 0)
    {
        YMIX::LogFile cf;
        cf << "Circuit " << name_ << ": Life time is " << scientific << setprecision(3) << timer_.get_dur() << " ms\n";
    }
}


void QCircuit::create_circ_file()
{
    cfname_ = path_to_output_ + "/" + name_ + FORMAT_CIRCUIT;
    if(env_.rank == 0)
    {
        YMIX::File cf(cfname_, true);
        cf << name_ << " " << nq_ << "\n";
    }
}

void QCircuit::create_tex_file()
{
    texname_ = path_to_output_ + "/" + name_ + FORMAT_TEX;
    if(env_.rank == 0)
    {
        YMIX::File cf(texname_, true);

        cf << "\\documentclass{article}\n";
        // cf << "\\usepackage[utf8]{inputenc}\n";
        cf << "\\usepackage[dvipsnames]{xcolor}\n";
        cf << "\\usepackage{amsmath}\n";
        cf << "\\usepackage{amsfonts,amssymb}\n";
        cf << "\\usepackage{tikz}\n";
        cf << "\\usetikzlibrary{quantikz}\n";
        cf << "\n";
        cf << "\\newcommand{\\yi}{\\mathrm{i}}\n";
        cf << "\n";
        cf << "\\begin{document}\n";
        cf << "\n\n";
        cf << "\\begin{figure}[t!]\n";
        cf << "\\centering\n";
        cf << "\\begin{quantikz}[row sep={0.5cm,between origins},column sep={0.3cm}]\n";   
    } 
}


void QCircuit::finish_tex_file()
{
    if(env_.rank == 0)
    {
        YMIX::File cf(texname_, false);

        // --- write the matrix of strings to the file ---
        int counter_row = -1;
        for(auto const & one_row_vec: tex_lines_)
        {
            counter_row++;
            string line_row = "";
            for(auto const& one_phrase: one_row_vec)
            {
                line_row += one_phrase;
            }   
            if(counter_row < tex_lines_.size())
                line_row += "\\\\"s;
            cf << line_row << "\n";
        }

        // --- finish the .tex file ---
        cf << "\\end{quantikz}\n";
        cf << "\\caption{Figure of the circuit: " << name_ << "}\n";
        cf << "\\end{figure}\n";
        cf << "\n\n";
        cf << "\\end{document}\n";
    }
}



void QCircuit::generate(const bool& flag_print)
{
    YMIX::YTimer timer;

    auto start = gates_.begin() + id_start_;
    for(auto it = start; it != gates_.end(); ++it)
    {
        // // // if(env_.rank == 0) cout << env_.rank << " Generating " << (*it)->get_name() << endl;
        // if(YMIX::compare_strings((*it)->get_name(), "cond_R") || YMIX::compare_strings((*it)->get_name(), "cond_R*"))
        //     timer.StartPrint(env_, "CondR generation... ");
        // if(YMIX::compare_strings((*it)->get_name(), "H"))
        //     timer.StartPrint(env_, "H generation... ");

        (*it)->generate(c_);

        // if(YMIX::compare_strings((*it)->get_name(), "cond_R") || YMIX::compare_strings((*it)->get_name(), "cond_R*"))
        //     timer.StopPrint(env_);
        // if(YMIX::compare_strings((*it)->get_name(), "H"))
        //     timer.StopPrint(env_);
    }
    id_start_ += gates_.end() - start;

    print_gates(flag_print); 
}

void QCircuit::generate(string& stop_name, int& id_current, const bool& flag_print)
{
    auto start  = gates_.begin() + id_start_;
    auto it_end = gates_.end();
    stop_name = name_;
    for(auto it = start; it != gates_.end(); ++it) 
    {
        if(YMIX::compare_strings((*it)->get_type(), "stop"))
        {
            it_end = it+1;
            stop_name = (*it)->get_name();
            break;
        }

        (*it)->generate(c_);
    }
    id_start_ += it_end - start;
    id_current = id_start_;

    print_gates(flag_print);
}

void QCircuit::conjugate_transpose()
{
    reverse(gates_.begin(), gates_.end());
    for(auto& gate: gates_)
    {
        //cout << "Inversing " << gate->get_name() << endl;
        gate->conjugate_transpose();
        gate->set_layer(
            oo_layers_->get_n_layers() - gate->get_layer()
        );
    }  
}

void QCircuit::print_gates(const bool& flag_print)
{
    if(env_.rank == 0)
    {
        if(flag_print)
        {
            // --- print gates to the .circuit file ---
            YMIX::File cf(cfname_);
            for(auto& gate: gates_)
                gate->write_to_file(cf);

            // // --- print gates to the .tex file ---
            // for(auto& gate: gates_)
            // {
            //     gate->get_tex_representation();


            // }
        }
    }
}

void QCircuit::copy_gates_from(YCCQ c, YCVI regs_new, YCCB box, YCB flag_inv)
{
    if(box)
    {
        YSG oo = box->copy_gate();
        oo_layers_->add_gate(oo);
        gates_.push_back(oo);
    }

    if(flag_inv)
    {
        auto gates_c = vector<YSG>(c->gates_);
        reverse(gates_c.begin(), gates_c.end());
        for(auto& gate: gates_c)
        {
            auto gate_copy = gate->copy_gate();
            gate_copy->conjugate_transpose();
            gate_copy->correct_qubits(regs_new);
            oo_layers_->add_gate(gate_copy);
            gates_.push_back(gate_copy);
        }
    }
    else
    {
        for(auto& gate: c->gates_)
        {
            auto gate_copy = gate->copy_gate();
            gate_copy->correct_qubits(regs_new);
            oo_layers_->add_gate(gate_copy);
            gates_.push_back(gate_copy);
        }
    }

    if(box)
    {
        YSG oo = box->copy_box();
        oo->set_flag_start(false);
        oo_layers_->add_gate(oo);
        gates_.push_back(oo);
    }
}

void QCircuit::insert_gates_from(const QCircuit* of, YCCB box)
{
    if(box)
    {
        YSG oo = box->copy_gate();
        oo_layers_->add_gate(oo);
        gates_.push_back(oo);
    }

    for(const auto& gate: of->gates_)
    {
        // !!! do not correct here the gate layers !!!
        gates_.push_back(gate);
    }
        
    if(box)
    {
        YSG oo = box->copy_box();
        oo->set_flag_start(false);
        oo_layers_->add_gate(oo);
        gates_.push_back(oo);
    }
}

YVIv QCircuit::add_register(YCS name, YCU n_qubits, YCB flag_ancilla)
{
    int nq_tot = c_.numQubitsRepresented;

    unsigned shift = 0;
    for(auto const& reg_name: regnames_)
        shift += regs_[reg_name].size();

    regnames_.push_back(name);

    // 0-th qubit is the least significant in the register;
    // the last qubit is the most significant:
    YVIv qubits_positions(n_qubits);
    for(unsigned i = 0; i < n_qubits; ++i)
        qubits_positions[i] = nq_tot - shift - n_qubits + i; // YINV

    regs_[name] = qubits_positions;

    if(flag_ancilla)
        ancs_.insert(
            ancs_.begin(), 
            qubits_positions.begin(), 
            qubits_positions.end()
        );

    return qubits_positions;
}

void QCircuit::print_reg_positions(std::ofstream& of) const
{
    for(auto const& reg_name: regnames_)
    {
        of << reg_name << " ";
        for(auto& id_q: regs_.at(reg_name))
            of << id_q << " ";
    }
}

void QCircuit::set_standart_output_format()
{
    int reg_size;
    int reg_size_total = 0;
    if(regs_.empty())
    {
        standart_output_format_.push_back(nq_);
        return;
    }
    
    for(auto const& reg_name: regnames_)
    {
        reg_size = regs_[reg_name].size();
        reg_size_total += reg_size;
        standart_output_format_.push_back(reg_size);
    }

    int diff = nq_ - reg_size_total;
    if(diff > 0) standart_output_format_.push_back(diff);
}

void QCircuit::save_regs()
{
    if(env_.rank == 0) 
    {
        if(regnames_.size() > 0)
        {
            YMIX::File cf(cfname_);
            cf << "QubitRegisterNames ";
            for(auto const& reg_name: regnames_)
                cf << reg_name << " " << regs_[reg_name].size() << " ";
            cf.of << endl;
        }

        for(auto const& [reg_name, reg_qs] : regs_)
        {
            auto reg_nq = reg_qs.size();
            for(auto i = 0; i < reg_nq; i++)
            {
                tex_lines_[i].push_back(
                    "\\lstick{$" + reg_name + "_" + to_string(reg_nq - i - 1) + "$}"
                );
            }
        }
    }
}

YQCP QCircuit::get_the_circuit()
{
    YQCP temp = this;
    return this;
}

void QCircuit::empty_binary_states()
{
    ib_state_ = vector<short>(nq_);
}
void QCircuit::prepare_zero_init_state()
{
    initZeroState(c_);
}
void QCircuit::reset()
{
    gates_.clear();
    id_start_ = 0;

    destroyQureg(c_, env_);
    c_ = createQureg(nq_, env_);

    if(!init_state_.empty())
        for(auto& state: init_state_)
            reset_init_vector(state);
    else{
        if(!ib_state_.empty())
            set_init_binary_state();
    }
}
void QCircuit::reset_qureg()
{
    id_start_ = 0;
    destroyQureg(c_, env_);
    c_ = createQureg(nq_, env_);
}

void QCircuit::set_init_binary_state(const bool& flag_mpi_bcast)
{
    if(flag_mpi_bcast)
    {
        int nq = c_.numQubitsRepresented;
        if(env_.rank > 0)
            if(ib_state_.empty())
                ib_state_ = vector<short>(nq);
        if(YMPI) MPI_Bcast(&ib_state_[0], nq, MPI_SHORT, 0, MPI_COMM_WORLD);
    }

    long long int ii = YMATH::binaryToInt(ib_state_);
    initClassicalState(c_, ii);
}

/**
 * !!! ATTENTION !!!: seem to be incorrect, qb odes not give correct position of an element within a state vector
 */
void QCircuit::set_init_vector(YCI qb, YCI nq, YVQ ampl_vec_real, YVQ ampl_vec_imag)
{
    long long b_ampl = 1 << qb;
    if(qb == 0) b_ampl = 0;

    long long n_ampls = 1 << nq;
    setAmps(c_, b_ampl, &ampl_vec_real[0], &ampl_vec_imag[0], n_ampls);

    INIT_STATE__ init_state;
    init_state.flag_defined = true;
    init_state.b_ampl = b_ampl;
    init_state.n_ampls = n_ampls;
    init_state.ampl_vec_real = YVQv(ampl_vec_real);
    init_state.ampl_vec_imag = YVQv(ampl_vec_imag);
    init_state_.push_back(init_state);
}

/**
 * !!! ATTENTION !!!: by adding init_state one by one, one creates a superposition.
 * Moreover, since at the beginning one has a zero state, this zero state is saved.
 * There is no renormalization in this method.
 */
void QCircuit::set_init_vector(INIT_STATE__& init_state)
{
    setAmps(
        c_, 
        init_state.b_ampl, 
        &init_state.ampl_vec_real[0], 
        &init_state.ampl_vec_imag[0], 
        init_state.n_ampls
    );
    init_state_.push_back(init_state);
}
void QCircuit::reset_init_vector(INIT_STATE__& state)
{
    setAmps(
        c_, state.b_ampl, 
        &state.ampl_vec_real[0], &state.ampl_vec_imag[0], 
        state.n_ampls
    );
}

void QCircuit::set_qubit_state(YCU id_q)
{
    if(ib_state_.empty())
        ib_state_ = vector<short>(nq_);
        
    if(id_q >= nq_)
    {
        cerr << "Error while setting init. state of qubits:\n";
        cerr << "A qubit with id = " << id_q << " is requested while there are only " 
            << nq_ << " qubits in a circuit " << name_ << "\n";
        cerr << "Remark: qubit indices start from zero." << endl;
        exit(-1);
    }
    ib_state_[nq_ - id_q - 1] = 1;

    // vector<short> bInt(nq_);
    // bInt[id_q] = 1;
    // auto b_ampl = YMATH::binaryToInt(bInt);

    // if(env_.rank == 0) cout << "resulting integer = " << b_ampl << endl;

    // INIT_STATE__ init_state;
    // init_state.flag_defined = true;
    // init_state.b_ampl = b_ampl;
    // init_state.n_ampls = 1;
    // init_state.ampl_vec_real = YVQv {1};
    // init_state.ampl_vec_imag = YVQv {0};
    // set_init_vector(init_state);
}

void QCircuit::set_qubit_state(YCVI ids_qs)
{
    if(ib_state_.empty())
        ib_state_ = vector<short>(nq_);

    for(auto const& id_q: ids_qs)
        set_qubit_state(id_q);
}

void QCircuit::set_reg_state(YCS name, YCI id_reg_qubit)
{
    if(ib_state_.empty())
        ib_state_ = vector<short>(nq_);

    if(YMIX::compare_strings(name, gv.reg_whole_circuit))
    {
        set_qubit_state(id_reg_qubit);
        return;
    }

    if(regs_.find(name) == regs_.end())
    {
        cerr << "\nError: No register with a name " << name << " in a circuit " << name_ << "." << endl;
        exit(-1);
    }

    auto reg_qubits = regs_[name];
    if(id_reg_qubit >= reg_qubits.size())
    {
        cerr << "\nError: The register (" << name << ") has only "
            << reg_qubits.size() << " qubits,\n";
        cerr << "while a qubit with id = " << id_reg_qubit  
            << " has been requested." << endl;
        exit(-1);
    }

    int id_circuit_qubit = reg_qubits[id_reg_qubit];

    // set_qubit_state(id_circuit_qubit);
    ib_state_[nq_ - id_circuit_qubit - 1] = 1;
}

void QCircuit::set_reg_state(YCS name, YCVI ids_reg_qubits)
{
    if(YMIX::compare_strings(name, gv.reg_whole_circuit))
    {
        set_qubit_state(ids_reg_qubits);
        return;
    }

    if(ib_state_.empty())
        ib_state_ = vector<short>(nq_);

    for(auto& id_reg_qubit: ids_reg_qubits)
        set_reg_state(name, id_reg_qubit);
}

YQCP QCircuit::x(YCVI ts, YVIv cs)
{ 
    for(auto& t:ts) x(t, cs);
    return get_the_circuit();
}
YQCP QCircuit::y(YCVI ts, YVIv cs)
{ 
    for(auto& t:ts) y(t, cs);
    return get_the_circuit();
}
YQCP QCircuit::z(YCVI ts, YVIv cs)
{ 
    for(auto& t:ts) z(t, cs);
    return get_the_circuit();
}
YQCP QCircuit::h(YCVI ts, YVIv cs)
{
    for(auto& t:ts) h(t, cs);
    return get_the_circuit();
}

void QCircuit::wavefunction_standard_analysis(
    const vector<unsigned>& organize_state, 
    const unsigned& prec
){
    if(env_.rank == 0) printf("\n--- Circuit %s: wavefunction standard analysis ---\n", name_.c_str());

    string str_wv;
    string str_wv_nz;
    get_wavefunction(organize_state, str_wv, str_wv_nz, prec);
    if(env_.rank == 0){
        cout << "Full wavefunction :\n";
        cout << str_wv << "\n";
        cout << "States with non-zero amplitudes:\n";
        cout << str_wv_nz << "\n";
    }
}

void QCircuit::get_wavefunction(
    const vector<unsigned>& organize_state, 
    YS str_wv, 
    YS str_wv_nz,
    const unsigned& prec
){
    vector<Complex> ampls;
    list<vector<short>> states;

    YMIX::Wavefunction(c_, str_wv, states, ampls, organize_state, prec);
        
    list<vector<short>> states_nz;
    vector<Complex> ampls_nz;
    YMIX::getNonzeroWavefunction(
        states, ampls, organize_state, 
        str_wv_nz, states_nz, ampls_nz,
        prec
    );
}

void QCircuit::get_special_wavefunction(
    const std::vector<short>& state_to_choose,
    const std::vector<unsigned>& organize_state, 
    YS str_wv_chosen,
    const unsigned& prec
){
    vector<Complex> ampls;
    list<vector<short>> states;
    string str_wv, str_wv_nz;
    YMIX::Wavefunction(c_, str_wv, states, ampls, organize_state, prec);
        
    list<vector<short>> states_nz;
    vector<Complex> ampls_nz;
    YMIX::getNonzeroWavefunction(
        states, ampls, organize_state, 
        str_wv_nz, states_nz, ampls_nz,
        prec
    );

    list<vector<short>> states_sp;
    vector<Complex> ampls_sp;
    YMIX::getSpecialStates(
        state_to_choose, states_nz, ampls_nz, 
        organize_state, 
        str_wv_chosen, 
        states_sp, ampls_sp, 
        prec 
    );
}
void QCircuit::get_special_wavefunction(
    const std::vector<short>& state_to_choose,
    const std::vector<unsigned>& organize_state, 
    YS str_wv_chosen,
    list<vector<short>>& states_chosen,
    vector<Complex>& ampls_chosen,
    const unsigned& prec
){
    vector<Complex> ampls;
    list<vector<short>> states;
    string str_wv, str_wv_nz;
    YMIX::Wavefunction(c_, str_wv, states, ampls, organize_state, prec);
        
    list<vector<short>> states_nz;
    vector<Complex> ampls_nz;
    YMIX::getNonzeroWavefunction(
        states, ampls, organize_state, 
        str_wv_nz, states_nz, ampls_nz,
        prec
    );

    YMIX::getSpecialStates(
        state_to_choose, states_nz, ampls_nz, 
        organize_state, 
        str_wv_chosen, states_chosen, ampls_chosen, 
        prec 
    );
}

void QCircuit::get_state_vector(YVQ state_real, YVQ state_imag)
{
    YMIX::get_state_vector(c_, nq_, state_real, state_imag);
}

void QCircuit::get_state_zero_ancillae(
    YCVU organize_state, 
    YS str_wv_out,
    list<vector<short>>& states_out,
    vector<Complex>& ampls_out,
    YCVsh state_to_choose,
    YCU prec
){
    unsigned n_states = nq_ - ancs_.size();
    YMIX::Wavefunction_NonzeroProbability(
        c_, n_states, organize_state, 
        str_wv_out, states_out, ampls_out, 
        state_to_choose, prec  
    );
}
void QCircuit::get_state_full(
    YCVU organize_state, 
    YS str_wv_out,
    list<vector<short>>& states_out,
    vector<Complex>& ampls_out,
    YCVsh state_to_choose,
    YCU prec
){
    unsigned n_states = nq_;
    YMIX::Wavefunction_NonzeroProbability(
        c_, n_states, organize_state, 
        str_wv_out, states_out, ampls_out, 
        state_to_choose, prec  
    );
}

void QCircuit::controlled(YCVI cs)
{
    for(auto& gate: gates_)
        gate->add_control_qubits(cs);
}
void QCircuit::controlled(YCI c)
{
    for(auto& gate: gates_)
        gate->add_control_qubits(c);
}

void QCircuit::read_structure_gate(
    YISS istr, YVI ids_target, qreal& par_gate, 
    YVI ids_control, YVI ids_x, 
    YVVI ids_control_it, YVVI ids_x_it
){
    string word;

    // --- read target qubits ---
    read_reg_int(istr, ids_target);
    
    // --- read the parameter of the gate ---     
    if(!isnan(par_gate))
    {
        istr >> word;
        par_gate = get_value_from_word(word);
    } 

    // --- read the end of the gate structure description ---
    read_end_gate(istr, ids_control, ids_x, ids_control_it, ids_x_it);  
}

void QCircuit::read_structure_gate(
    YISS istr, YVI ids_target, qreal& par_gate1, qreal& par_gate2,
    YVI ids_control, YVI ids_x, 
    YVVI ids_control_it, YVVI ids_x_it
){
    string word;

    // --- read target qubits ---
    read_reg_int(istr, ids_target);
    
    // --- read parameters of the gate ---     
    if(!isnan(par_gate1))
    {
        istr >> word;
        par_gate1 = get_value_from_word(word);
    } 
    if(!isnan(par_gate2))
    {
        istr >> word;
        par_gate2 = get_value_from_word(word);
    } 

    // --- read the end of the gate structure description ---
    read_end_gate(istr, ids_control, ids_x, ids_control_it, ids_x_it);  
}

void QCircuit::read_end_gate(YISS istr, YVI ids_control, YVI ids_x, YVVI ids_control_it, YVVI ids_x_it)
{
    std::string word;
    while(istr >> word)
    {
        if(YMIX::compare_strings(word, "end_gate"))
            return;

        if(YMIX::compare_strings(word, "control"))
            read_reg_int(istr, ids_control);

        if(YMIX::compare_strings(word, "ocontrol"))
        {
            YVIv ids_c_local;
            read_reg_int(istr, ids_c_local);
            ids_x.insert(
                ids_x.end(), 
                ids_c_local.begin(), 
                ids_c_local.end()
            );
            ids_control.insert(
                ids_control.end(), 
                ids_c_local.begin(), 
                ids_c_local.end()
            );
        }

        if(YMIX::compare_strings(word, "control_it"))
        {
            YVIv ids_local;
            read_reg_int(istr, ids_local);
            ids_control_it.push_back(ids_local);
        } 

        if(YMIX::compare_strings(word, "ocontrol_it"))
        {
            YVIv ids_local;
            read_reg_int(istr, ids_local);
            ids_control_it.push_back(ids_local);
            ids_x_it.push_back(ids_local);
        }
    }
}

void QCircuit::read_reg_int(YISS istr, YVI ids_target, YCS word_start)
{
    string reg_name, word;
    bool flag_read_reg_name = false;
    int n_regs, integer_qu, n_bitA;
    int nq_reg;

    // nodes can sit on several registers
    if(word_start.empty())
        istr >> word;
    else 
        word = word_start;

    if(YMATH::is_number(word))
    {
        n_regs = stoi(word);
        flag_read_reg_name = true;
    }
    else
    {
        n_regs = 1;
        reg_name = word;
        flag_read_reg_name = false;
    }

    // within every register, one can have several qubits
    for(unsigned i_reg = 0; i_reg < n_regs; i_reg++)
    {
        if(flag_read_reg_name) istr >> reg_name;
        if(!YMIX::is_present(regnames_, reg_name))
        {
            if(YMPI) MPI_Barrier(MPI_COMM_WORLD);
            throw "no register with the name " + reg_name;
            if(YMPI) MPI_Barrier(MPI_COMM_WORLD);
        }

        try
        {
            istr >> word;
            integer_qu = get_value_from_word(word);

            auto reg_chosen = regs_[reg_name];
            nq_reg = reg_chosen.size();
            if(integer_qu < 0)
                integer_qu = (1 << nq_reg) + integer_qu;

            vector<short> binArray(nq_reg);
            YMATH::intToBinary(integer_qu, binArray);

            n_bitA = binArray.size();
            for(unsigned id_bit = 0; id_bit < n_bitA; id_bit++)
                if(binArray[n_bitA - id_bit - 1] == 1)
                    ids_target.push_back(reg_chosen[id_bit]);
        }
        catch(YCS e)
        {
            throw "wrong format of the number of qubits in the register " + reg_name + ": " + e;
        }
    }
}

void QCircuit::read_structure_gate_condR_split(YISS istr, YCS path_in, YCB flag_inv)
{
    string name;
    YVIv ids_target, ids_cond, ids_control, ids_x;
    int int_sel_flag = 0;
    long long size_cond;
    vector<qreal> as; // profile of rotation angles;

    // --- name of the conditional rotation gate ---
    istr >> name;

    // --- read target qubits ---
    read_reg_int(istr, ids_target);

    // --- read conditional qubits ---
    read_reg_int(istr, ids_cond);
    size_cond = pow(2, ids_cond.size());
    as = YVQv(size_cond);

    // --- read end of gate structure ---
    YVVIv ids_control_it, ids_x_it;
    read_end_gate(istr, ids_control, ids_x, ids_control_it, ids_x_it);

    // -------------------------------------------
    // --- read a file with condR parameters ---
    string data;
    string ifname = path_in + "/" + name + FORMAT_PROFILE;
    char cdata[YSIZE_CHAR_ARRAY];
    if(env_.rank == 0)
    {
        ifstream ff(ifname);
        if(!ff.is_open()) throw "Error: there is not a file: " + ifname;
        data = string((istreambuf_iterator<char>(ff)), istreambuf_iterator<char>());
        ff.close();
    }

    int size_data;
    if(env_.rank == 0) 
    {
        if(YSIZE_CHAR_ARRAY < data.size()) 
            throw "Error: Size of the char array is too small to save the input file."s;

        strcpy(cdata, data.c_str());
        size_data = data.size();
    }
    if(YMPI) MPI_Bcast(&size_data, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
    if(YMPI) MPI_Bcast(cdata, size_data+1, MPI_CHAR, 0, MPI_COMM_WORLD);
    if(env_.rank > 0) data = string(cdata);

    istringstream iss(data);

    // --- read value profiles ---
    qreal vv;
    for(unsigned i = 0; i < size_cond; ++i)
    {
        iss >> vv; // value from a profile that has to be normalized to 1
        as[i] = 2.0 * acos(vv); // angle of rotation
    }

    // -----------------------------------------
    // --- create conditional rotation gates ---
    sort(ids_cond.begin(),    ids_cond.end());
    sort(ids_control.begin(), ids_control.end());

    // add 0-control
    x(ids_x);

    // add a set of Ry gates controlled on conditional (and control) qubits:
    unsigned N_cond = 1 << ids_cond.size();
    unsigned count_mod;
    vector<int> full_c = vector<int>(ids_control);
    full_c.insert(full_c.end(), ids_cond.begin(), ids_cond.end());

    for(auto const& id_target: ids_target) 
    {
        x(ids_cond);
        ry(id_target, as[0], full_c, flag_inv);

        for(unsigned id_as = 1; id_as < N_cond; id_as++)
        {
            x(ids_cond[0]);
            count_mod = 1;
            for(unsigned int_mod = 2; int_mod < N_cond; int_mod *= 2)
            {
                if(id_as % int_mod == 0) x(ids_cond[count_mod]);
                count_mod++;
            }
            ry(id_target, as[id_as], full_c, flag_inv);
        }
    }
        
    // close 0-control
    x(ids_x);
}

void QCircuit::read_structure_gate_adder1(YISS istr, YCS path_in, YCB flag_inv)
{
    YVIv ids_target, ids_control, ids_x;
    long long nt;

    // --- read target qubits ---
    read_reg_int(istr, ids_target);

    // --- read end of gate structure ---
    YVVIv ids_control_it, ids_x_it;
    read_end_gate(istr, ids_control, ids_x, ids_control_it, ids_x_it);

    // --- put the high-priority qubits at the beginning ---
    sort(ids_target.begin(), ids_target.end(), greater<int>());
    nt = ids_target.size();

    // --- add CNOT and X gates with control nodes ---
    x(ids_x);

    for(unsigned i = 0; i < nt-1; ++i)
    {
        YVIv ids_cnot_cs = YVIv(ids_target.begin() + i + 1, ids_target.end());
        ids_cnot_cs.insert(ids_cnot_cs.end(), ids_control.begin(), ids_control.end());
        x(ids_target[i], ids_cnot_cs);
    }
    x(ids_target.back(), ids_control);

    x(ids_x);
}

void QCircuit::read_structure_gate_subtractor1(YISS istr, YCS path_in, YCB flag_inv)
{
    YVIv ids_target, ids_control, ids_x;
    long long nt;

    // --- read target qubits ---
    read_reg_int(istr, ids_target);

    // --- read end of gate structure ---
    YVVIv ids_control_it, ids_x_it;
    read_end_gate(istr, ids_control, ids_x, ids_control_it, ids_x_it);

    // --- put the low-priority qubits at the beginning ---
    sort(ids_target.begin(), ids_target.end());
    nt = ids_target.size();

    // --- add CNOT and X gates with control nodes ---
    x(ids_x);
    x(ids_target[0], ids_control);
    for(unsigned i = 1; i < nt; ++i)
    {
        YVIv ids_cnot_cs = YVIv(ids_target.begin(), ids_target.begin() + i);
        ids_cnot_cs.insert(ids_cnot_cs.end(), ids_control.begin(), ids_control.end());
        x(ids_target[i], ids_cnot_cs);
    }
    x(ids_x);
}

void QCircuit::read_structure_gate_swap(YISS istr, YCS path_in, YCB flag_inv)
{
    YVIv ids_target_1, ids_target_2, ids_control, ids_x;
    long long nt;

    // --- read target qubits ---
    read_reg_int(istr, ids_target_1);
    read_reg_int(istr, ids_target_2);

    // --- read end of gate structure ---
    YVVIv ids_control_it, ids_x_it;
    read_end_gate(istr, ids_control, ids_x, ids_control_it, ids_x_it);

    // --- add CNOT gates ---
    nt = ids_target_1.size();

    x(ids_x);
    for(unsigned i = 0; i < nt; ++i)
        swap(ids_target_1[i], ids_target_2[i], ids_control);
    x(ids_x);
}

YQCP QCircuit::swap(YCI t1, YCI t2, YVIv cs)
{ 
    YVIv ids_cs_1 = {t1};
    YVIv ids_cs_2 = {t2};
    ids_cs_1.insert(ids_cs_1.end(), cs.begin(), cs.end());
    ids_cs_2.insert(ids_cs_2.end(), cs.begin(), cs.end());
    return x(t2, ids_cs_1)->x(t1, ids_cs_2)->x(t2, ids_cs_1); 
}

qreal QCircuit::get_value_from_word(YCS word)
{
    if(word.find("<") == string::npos)
    {
        istringstream sstr(word);
        qreal res_value;
        if(!(sstr >> res_value))
            throw "Wrong format"s;
        return res_value;
    }

    unsigned first = word.find("<");
    unsigned last = word.find(">");
    string const_name = word.substr(first+1,last-first-1);

    if(constants_.find(const_name) == constants_.end())
        throw "A constant with the name "s + const_name + " is not found."s;
    return constants_[const_name];
}