#include "../include/QCircuit.h"

using namespace std;




QCircuit::QCircuit(
    YCS name, 
    const QuESTEnv& env, 
    YCS path_to_output, 
    YCU nq, 
    const std::map<std::string, qreal>& constants,
    YCB flag_circuit,
    YCB flag_tex,
    YCB flag_layers
) :
name_(name),
env_(env),
path_to_output_(path_to_output),
flag_circuit_(flag_circuit),
flag_tex_(flag_tex),
flag_layers_(flag_layers)
{
    timer_.Start();
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

    tex_noc_.resize(nq_, 1);
    tex_lines_.resize(nq_);

    create_circ_file();
    create_tex_file();
    initZeroState(c_);

    if(flag_layers_) oo_layers_ = make_shared<CircuitLayers__>(nq_);
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
    init_state_     = oc->init_state_;
    regs_           = map<string, YVIv>(oc->regs_);
    flags_anc_regs_ = map<string, bool>(oc->flags_anc_regs_);
    regnames_       = vector<string>(oc->regnames_);
    ancs_           = YVIv(oc->ancs_);
    ib_state_       = vector<short>(oc->ib_state_);
    id_start_       = oc->id_start_;
    standart_output_format_ = oc->standart_output_format_;
    unique_gates_names_ = YVSv(oc->unique_gates_names_);

    flag_circuit_ = oc->flag_circuit_;
    flag_tex_     = oc->flag_tex_;
    flag_layers_  = oc->flag_layers_;

    tex_noc_ = vector<uint64_t>(oc->tex_noc_);
    tex_lines_ = vector<vector<string>>(oc->tex_lines_);

    create_circ_file();
    create_tex_file();
    save_regs();

    for(auto const& gate: oc->gates_)
    {
        auto gate_copy = gate->copy_gate();
        gates_.push_back(gate_copy);
    }
    constants_ = map<string, qreal>(oc->constants_);
    if(flag_layers_) oo_layers_ = make_shared<CircuitLayers__>(oc->oo_layers_);
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
    if(!flag_circuit_)
        return;

    cfname_ = path_to_output_ + "/" + name_ + FORMAT_CIRCUIT;
    if(env_.rank == 0)
    {
        YMIX::File cf(cfname_, true);
        cf << name_ << " " << nq_ << "\n";
    }
}


void QCircuit::create_tex_file()
{
    if(!flag_tex_)
    {
        // cout << "HERE: circuit [" << name_ << "]: no tex" << endl;
        return;
    }
        

    texname_ = path_to_output_ + "/" + name_ + FORMAT_TEX;
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
    cf << TEX_BEGIN_CIRCUIT;   
}


void QCircuit::finish_tex_file()
{
    if(!flag_tex_)
        return;

    if(env_.rank == 0)
    {
        YMIX::File cf(texname_, false);

        // --- adjust the widths of the .tex phrases ---
        auto n_layers = tex_lines_[0].size();
        uint32_t width_max, width_curr;
        vector<uint32_t> q_widths(nq_);
        for(auto id_layer = 1; id_layer < n_layers; id_layer++)
        {
            if(id_layer%YGV::tex_circuit_length == 0)
            {
                for(auto id_q = 0; id_q < nq_; id_q++)
                    q_widths[id_q] = 0;
                continue;
            }

            for(auto id_q = 0; id_q < nq_; id_q++)
                q_widths[id_q] += tex_lines_[id_q][id_layer-1].length();

            width_max = 0;
            for(auto id_q = 0; id_q < nq_; id_q++)
                if(width_max < q_widths[id_q])
                    width_max = q_widths[id_q];

            for(auto id_q = 0; id_q < nq_; id_q++)
            {
                width_curr = q_widths[id_q];
                string line_padding(width_max - width_curr, ' ');

                tex_lines_[id_q][id_layer] = line_padding + tex_lines_[id_q][id_layer];
            }
        }

        // --- write the matrix of strings to the file ---
        int n_circ_pieces = int(tex_lines_[0].size() / YGV::tex_circuit_length)+1;
        if(tex_lines_[0].size() % YGV::tex_circuit_length == 0)
            n_circ_pieces--;
        for(auto id_piece = 0; id_piece < n_circ_pieces; id_piece++)
        {
            int counter_row = -1;
            int column_begin = id_piece * YGV::tex_circuit_length;
            int column_end = std::min(column_begin + YGV::tex_circuit_length, int(tex_lines_[0].size()));
            for(auto const & one_row_vec: tex_lines_)
            {
                counter_row++;
                string line_row = "";
                for(auto id_c = column_begin; id_c < column_end; id_c++)
                    line_row += one_row_vec[id_c];
                if(counter_row < (tex_lines_.size()-1))
                    line_row += "\\\\"s;
                cf << line_row << "\n";
            }
            cf << "\\end{quantikz}";

            if(id_piece < (n_circ_pieces - 1))
            {
                cf << "\\\\\n";
                cf << TEX_VERTICAL_GAP << "\n";
                cf << TEX_BEGIN_CIRCUIT;
            }
            else
            {
                cf << "\n";
            }
        }

        // --- finish the .tex file ---
        cf << "\\caption{Figure of the circuit: " << name_ << "}\n";
        cf << "\\end{figure}\n";
        cf << "\n\n";
        cf << "\\end{document}\n";
    }
}


void QCircuit::generate()
{
    YMIX::YTimer timer;

    auto start = gates_.begin() + id_start_;
    for(auto it = start; it != gates_.end(); ++it)
    {
        (*it)->generate(c_);
    }
    id_start_ += gates_.end() - start;
}


void QCircuit::generate(string& stop_name, int& id_current)
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
}


void QCircuit::print_gates()
{
    // -------------------------------------------------
    // --- print gates to the .circuit file ---
    if(env_.rank == 0 && flag_circuit_)
    {
        
        YMIX::File cf(cfname_);
        for(auto& gate: gates_)
            gate->write_to_file(cf);
    }

    // -------------------------------------------------
    // --- print gates to the .tex file ---
    if(env_.rank == 0 && flag_tex_)
    {
        YVIv ids_qubits_of_gate, ids_range, ids_targets, ids_controls;
        uint64_t id_first_noc_layer;
        uint64_t id_new_noc_layer;
        int id_b, id_t;
        bool flag_box = false;
        string box_name;
        for(auto& gate: gates_)
        {
            id_first_noc_layer = 0;
            if(YMIX::compare_strings(gate->get_type(), "stop"))
                continue;

            if(YMIX::compare_strings(gate->get_type(), "box"))
            {
                string box_name_curr = gate->get_name();
                if(gate->get_flag_start())
                {
                    if(!flag_box)
                        box_name = box_name_curr;
                    flag_box = true;
                }
                else
                {
                    if(YMIX::compare_strings(box_name, box_name_curr))
                        flag_box = false;
                }  
            }
                
            if(flag_box)
                continue;

            gate->get_gubits_act_on(ids_qubits_of_gate);
            id_t = *(max_element(ids_qubits_of_gate.begin(), ids_qubits_of_gate.end()));
            id_b = *(min_element(ids_qubits_of_gate.begin(), ids_qubits_of_gate.end()));

            // --- find the first layer in .tex, where there is enough free qubits to place the gate ---
            for(auto id_qubit = id_b; id_qubit <= id_t; id_qubit++)
            {
                if(id_first_noc_layer < tex_noc_[nq_ - id_qubit - 1])
                    id_first_noc_layer = tex_noc_[nq_ - id_qubit - 1];
            }

            // --- put the gate to the .tex layer ---
            gate->write_tex(tex_lines_, id_first_noc_layer, nq_);


            // cout << "\n\n";
            // cout << "id_first_noc_layer: " << id_first_noc_layer << "\n";
            // cout << "[" << id_b << ", " << id_t << "]\n";
            // int counter_row = -1;
            // for(auto const & one_row_vec: tex_lines_)
            // {
            //     counter_row++;
            //     string line_row = "";
            //     for(auto const& one_phrase: one_row_vec)
            //         line_row += one_phrase;
            //     cout << line_row << "\n";
            // }
            // cout << "--- current tex_noc ---\n";
            // for(auto id_q = 0; id_q < nq_; id_q++)
            //     cout << "tex_noc_[" << nq_ - id_q - 1 << "] = " << tex_noc_[id_q] << "\n";


            // --- shift the ids of non-occupied layers in the .tex ---
            id_new_noc_layer = id_first_noc_layer + 1;
            for(auto id_qubit = id_b; id_qubit <= id_t; id_qubit++)
                tex_noc_[nq_ - id_qubit - 1] = id_new_noc_layer;

            // --- if necessary, add next empty .tex layer ---
            if(id_new_noc_layer >= tex_lines_[0].size())
                for(auto id_q = 0; id_q < nq_; id_q++)
                    tex_lines_[id_q].push_back("&\\qw"s);
        }
    }
}


void QCircuit::conjugate_transpose()
{
    reverse(gates_.begin(), gates_.end());
    for(auto& gate: gates_)
    {
        gate->conjugate_transpose();
        if(flag_layers_) gate->set_layer(oo_layers_->get_n_layers() - gate->get_layer());
    }  
}


void QCircuit::copy_gates_from(YCCQ c, YCVI regs_new, YCCB box, YCB flag_inv, YCVI cs)
{
    if(box)
    {
        YSG oo = box->copy_gate();
        if(!cs.empty()) oo->add_control_qubits(cs);
        if(flag_layers_) oo_layers_->add_gate(oo);
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
            if(!cs.empty()) gate_copy->add_control_qubits(cs);
            if(flag_layers_) oo_layers_->add_gate(gate_copy);
            gates_.push_back(gate_copy);
        }
    }
    else
    {
        for(auto& gate: c->gates_)
        {
            auto gate_copy = gate->copy_gate();
            gate_copy->correct_qubits(regs_new);
            if(!cs.empty()) gate_copy->add_control_qubits(cs);
            if(flag_layers_) oo_layers_->add_gate(gate_copy);
            gates_.push_back(gate_copy);
        }
    }

    if(box)
    {
        YSG oo = box->copy_box();
        oo->set_flag_start(false);
        if(!cs.empty()) oo->add_control_qubits(cs);
        if(flag_layers_) oo_layers_->add_gate(oo);
        gates_.push_back(oo);
    }

    unique_gates_names_.insert(
        unique_gates_names_.end(), 
        c->unique_gates_names_.begin(), 
        c->unique_gates_names_.end()
    );
}


void QCircuit::insert_gates_from(const QCircuit* c, YCCB box)
{
    if(box)
    {
        YSG oo = box->copy_gate();
        if(flag_layers_) oo_layers_->add_gate(oo);
        gates_.push_back(oo);
    }

    for(const auto& gate: c->gates_)
    {
        // !!! here, do not add the inserted gates to the layers !!!
        gates_.push_back(gate);
    }
        
    if(box)
    {
        YSG oo = box->copy_box();
        oo->set_flag_start(false);
        if(flag_layers_) oo_layers_->add_gate(oo);
        gates_.push_back(oo);
    }

    unique_gates_names_.insert(
        unique_gates_names_.end(), 
        c->unique_gates_names_.begin(), 
        c->unique_gates_names_.end()
    );
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
    {
        ancs_.insert(
            ancs_.begin(), 
            qubits_positions.begin(), 
            qubits_positions.end()
        );
        flags_anc_regs_[name] = true;
    }
    else
    {
        flags_anc_regs_[name] = false;
    }
        
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
    if(env_.rank == 0 && flag_circuit_) 
    {
        if(regnames_.size() > 0)
        {
            YMIX::File cf(cfname_);
            cf << "QubitRegisterNames ";
            for(auto const& reg_name: regnames_)
            {
                cf << reg_name << " " << regs_[reg_name].size() << " ";
                if(flags_anc_regs_[reg_name])
                    cf << "1" << " ";
                else
                    cf << "0" << " ";
            }
            cf.of << endl;
        }
    }

    if(env_.rank == 0 && flag_tex_) 
    {
        int counter_q = -1;
        for(auto const& reg_name: regnames_)
        // for(auto const& [reg_name, reg_qs] : regs_)
        {
            auto reg_qs = regs_[reg_name];
            auto reg_nq = reg_qs.size();

            string reg_name_part1, reg_name_part2;
            std::size_t pos = reg_name.find("_");
            if(pos != std::string::npos)
            {
                reg_name_part1 = reg_name.substr(0, pos);
                reg_name_part2 = reg_name.substr(pos+1) + ", ";
            }
            else
            {
                reg_name_part1 = string(reg_name);
                reg_name_part2 = "";
            }

            for(auto i = 0; i < reg_nq; i++)
            {
                counter_q++;

                // register name and qubit id within the register:
                tex_lines_[counter_q].push_back(
                    "\\lstick{$" + reg_name_part1 + "_{" + reg_name_part2 + to_string(reg_nq - i - 1) + "}$}"
                );

                // first layer:
                tex_lines_[counter_q].push_back("&\\qw");
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

    if(init_state_.flag_defined)
        reset_init_vector(init_state_);
    else{
        if(!ib_state_.empty())
            set_init_binary_state();
    }
    unique_gates_names_.clear();
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
    }

    long long int ii = YMATH::binaryToInt(ib_state_);
    initClassicalState(c_, ii);
}


void QCircuit::set_init_vector(YVQ ampl_vec_real, YVQ ampl_vec_imag)
{
    long long b_ampl = 0;
    long long N_ampls = ampl_vec_real.size();
    if(N_ampls != ampl_vec_imag.size())
        throw string("vectors with real and imaginary parts of the initial state are of different size.");
    setAmps(c_, b_ampl, &ampl_vec_real[0], &ampl_vec_imag[0], N_ampls);

    init_state_.flag_defined = true;
    init_state_.b_ampl = b_ampl;
    init_state_.n_ampls = N_ampls;
    init_state_.ampl_vec_real = YVQv(ampl_vec_real);
    init_state_.ampl_vec_imag = YVQv(ampl_vec_imag);
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

    if(YMIX::compare_strings(name, YGV::reg_whole_circuit))
    {
        set_qubit_state(id_reg_qubit);
        return;
    }

    if(regs_.find(name) == regs_.end())
    {
        cerr << "\nError: No register with the name " << name << " in a circuit " << name_ << "." << endl;
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
    ib_state_[nq_ - id_circuit_qubit - 1] = 1;
}

void QCircuit::set_reg_state(YCS name, YCVI ids_reg_qubits)
{
    if(YMIX::compare_strings(name, YGV::reg_whole_circuit))
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


void QCircuit::get_ref_to_state_vector(qreal*& state_real, qreal*& state_imag)
{
    copyStateFromGPU(c_);
    state_real = c_.stateVec.real;
    state_imag = c_.stateVec.imag;
}


void QCircuit::get_state(YMIX::StateVectorOut& out, YCB flag_ZeroPriorAnc)
{
    out.n_low_prior_qubits = flag_ZeroPriorAnc ? (nq_ - ancs_.size()): nq_;
    out.organize_state = get_standart_output_format();
    YMIX::Wavefunction_NonzeroProbability(c_, out);
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
    try
    {
        if(!isnan(par_gate))
        {
            istr >> word;
            par_gate = get_value_from_word(word);
        } 
    }   
    catch(YCS e)
    {
        throw "error in the format of the gate parameter: oracletool sees [" + word + "]: " + e;
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
    try
    {   
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
    }
    catch(YCS e)
    {
        throw "error in the format of the gate parameter: oracletool sees [" + word + "]: " + e;
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


void QCircuit::read_reg_int(YISS istr, YVI ids_target, YCB flag_sort, YCS word_start)
{
    string reg_name, word;
    bool flag_read_reg_name = false;
    int n_regs, integer_qu, n_bitA;
    int nq_reg;

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
            throw "no register with the name " + reg_name;

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

    // sort the final array:
    if(flag_sort)
        sort(ids_target.begin(), ids_target.end());
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

    ifstream ff(ifname);
    if(!ff.is_open()) throw "Error: there is not a file: " + ifname;
    data = string((istreambuf_iterator<char>(ff)), istreambuf_iterator<char>());
    ff.close();
    
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

    // --- read end of the gate structure ---
    YVVIv ids_control_it, ids_x_it;
    read_end_gate(istr, ids_control, ids_x, ids_control_it, ids_x_it);

    // add the adder:
    x(ids_x);
    adder_1(ids_target, ids_control, flag_inv);
    x(ids_x);
}


void QCircuit::read_structure_gate_adder2(YISS istr, YCS path_in, YCB flag_inv)
{
    YVIv ids_target, ids_control, ids_x;
    long long nt;

    // --- read target qubits ---
    read_reg_int(istr, ids_target);

    // --- read end of the gate structure ---
    YVVIv ids_control_it, ids_x_it;
    read_end_gate(istr, ids_control, ids_x, ids_control_it, ids_x_it);

    // add the adder:
    x(ids_x);
    adder_2(ids_target, ids_control, flag_inv);
    x(ids_x);
}


void QCircuit::read_structure_gate_adder3(YISS istr, YCS path_in, YCB flag_inv)
{
    YVIv ids_target, ids_control, ids_x;
    long long nt;

    // --- read target qubits ---
    read_reg_int(istr, ids_target);

    // --- read end of the gate structure ---
    YVVIv ids_control_it, ids_x_it;
    read_end_gate(istr, ids_control, ids_x, ids_control_it, ids_x_it);

    // add the adder:
    x(ids_x);
    adder_3(ids_target, ids_control, flag_inv);
    x(ids_x);
}


void QCircuit::read_structure_gate_subtractor1(YISS istr, YCS path_in, YCB flag_inv)
{
    YVIv ids_target, ids_control, ids_x;

    // --- read target qubits ---
    read_reg_int(istr, ids_target);

    // --- read end of gate structure ---
    YVVIv ids_control_it, ids_x_it;
    read_end_gate(istr, ids_control, ids_x, ids_control_it, ids_x_it);

    // add the subtractor:
    x(ids_x);
    subtractor_1(ids_target, ids_control, flag_inv);
    x(ids_x);
}


void QCircuit::read_structure_gate_subtractor2(YISS istr, YCS path_in, YCB flag_inv)
{
    YVIv ids_target, ids_control, ids_x;

    // --- read target qubits ---
    read_reg_int(istr, ids_target);

    // --- read end of gate structure ---
    YVVIv ids_control_it, ids_x_it;
    read_end_gate(istr, ids_control, ids_x, ids_control_it, ids_x_it);

    // add the subtractor:
    x(ids_x);
    subtractor_2(ids_target, ids_control, flag_inv);
    x(ids_x);
}


void QCircuit::read_structure_gate_subtractor3(YISS istr, YCS path_in, YCB flag_inv)
{
    YVIv ids_target, ids_control, ids_x;

    // --- read target qubits ---
    read_reg_int(istr, ids_target);

    // --- read end of gate structure ---
    YVVIv ids_control_it, ids_x_it;
    read_end_gate(istr, ids_control, ids_x, ids_control_it, ids_x_it);

    // add the subtractor:
    x(ids_x);
    subtractor_3(ids_target, ids_control, flag_inv);
    x(ids_x);
}


void QCircuit::read_structure_gate_adder(YISS istr, YCS path_in, YCB flag_inv)
{
    YVIv ids_target_1, ids_target_2, ids_target_3, ids_control, ids_x;
    long long nt;

    // --- read target qubits ---
    read_reg_int(istr, ids_target_1);
    read_reg_int(istr, ids_target_2);
    read_reg_int(istr, ids_target_3);

    if(ids_target_1.size() != ids_target_2.size())
        throw string("target registers must have the same size");
    if(ids_target_1.size() != ids_target_3.size())
        throw string("target registers must have the same size");
    if(ids_target_2.size() != ids_target_3.size())
        throw string("target registers must have the same size");

    // --- read the end of the gate structure ---
    YVVIv ids_control_it, ids_x_it;
    read_end_gate(istr, ids_control, ids_x, ids_control_it, ids_x_it);

    // --- adder of two variables ---
    x(ids_x); 
    adder(ids_target_1, ids_target_2, ids_target_3, ids_control, flag_inv);
    x(ids_x);
}


void QCircuit::read_structure_gate_subtractor(YISS istr, YCS path_in, YCB flag_inv)
{
    YVIv ids_target_1, ids_target_2, ids_target_3, ids_control, ids_x;
    long long nt;

    // --- read target qubits ---
    read_reg_int(istr, ids_target_1);
    read_reg_int(istr, ids_target_2);
    read_reg_int(istr, ids_target_3);

    if(ids_target_1.size() != ids_target_2.size())
        throw string("target registers must have the same size");
    if(ids_target_1.size() != ids_target_3.size())
        throw string("target registers must have the same size");
    if(ids_target_2.size() != ids_target_3.size())
        throw string("target registers must have the same size");

    // --- read the end of the gate structure ---
    YVVIv ids_control_it, ids_x_it;
    read_end_gate(istr, ids_control, ids_x, ids_control_it, ids_x_it);

    // --- adder of two variables ---
    x(ids_x); 
    subtractor(ids_target_1, ids_target_2, ids_target_3, ids_control, flag_inv);
    x(ids_x);
}


void QCircuit::read_structure_gate_adder_qft(YISS istr, YCS path_in, YCB flag_inv)
{
    YVIv qs_v1, qs_v2, qs_carry, ids_control, ids_x;
    int q_carry;
    YVVIv ids_control_it, ids_x_it;
    long long nt;

    // --- read target qubits ---
    read_reg_int(istr, qs_v1);
    read_reg_int(istr, qs_v2);
    read_reg_int(istr, qs_carry);
    q_carry = qs_carry[0];

    if(qs_v1.size() != qs_v2.size())
        throw string("target registers must have the same size");

    // --- read the end of the gate structure ---
    read_end_gate(istr, ids_control, ids_x, ids_control_it, ids_x_it);

    // --- adder ---
    x(ids_x); 
    adder_qft(qs_v1, qs_v2, q_carry, ids_control, flag_inv);
    x(ids_x);
}


void QCircuit::read_structure_gate_subtractor_qft(YISS istr, YCS path_in, YCB flag_inv)
{
    YVIv qs_v1, qs_v2, qs_sign, ids_control, ids_x;
    int q_sign;
    YVVIv ids_control_it, ids_x_it;
    long long nt;

    // --- read target qubits ---
    read_reg_int(istr, qs_v1);
    read_reg_int(istr, qs_v2);
    read_reg_int(istr, qs_sign);
    q_sign = qs_sign[0];

    if(qs_v1.size() != qs_v2.size())
        throw string("target registers must have the same size");

    // --- read the end of the gate structure ---
    read_end_gate(istr, ids_control, ids_x, ids_control_it, ids_x_it);

    // --- adder ---
    x(ids_x); 
    subtractor_qft(qs_v1, qs_v2, q_sign, ids_control, flag_inv);
    x(ids_x);
}



void QCircuit::read_structure_gate_adder_fixed(YISS istr, YCS path_in, YCB flag_inv)
{
    YVIv ids_target, ids_carry_temp, ids_control, ids_x;
    int  id_carry;
    uint32_t int_sub;
    string word;

    // --- read the gate parameters ---
    read_reg_int(istr, ids_target);
    
    istr >> word;
    int_sub = get_value_from_word(word);

    read_reg_int(istr, ids_carry_temp);
    id_carry = ids_carry_temp[0];

    // --- read the end of the gate structure ---
    YVVIv ids_control_it, ids_x_it;
    read_end_gate(istr, ids_control, ids_x, ids_control_it, ids_x_it);

    // --- adder of two variables ---
    x(ids_x); 
    adder_fixed(ids_target, id_carry, int_sub, ids_control, flag_inv);
    x(ids_x);
}


void QCircuit::read_structure_gate_subtractor_fixed(YISS istr, YCS path_in, YCB flag_inv)
{
    YVIv ids_target, ids_carry_temp, ids_control, ids_x;
    int  id_carry;
    uint32_t int_sub;
    string word;

    // --- read the gate parameters ---
    read_reg_int(istr, ids_target);
    
    istr >> word;
    int_sub = get_value_from_word(word);

    read_reg_int(istr, ids_carry_temp);
    id_carry = ids_carry_temp[0];

    // --- read the end of the gate structure ---
    YVVIv ids_control_it, ids_x_it;
    read_end_gate(istr, ids_control, ids_x, ids_control_it, ids_x_it);

    // --- gate---
    x(ids_x); 
    subtractor_fixed(ids_target, id_carry, int_sub, ids_control, flag_inv);
    x(ids_x);
}


void QCircuit::read_structure_gate_comparator_fixed(YISS istr, YCS path_in, YCB flag_inv)
{
    YVIv ids_target, ids_carry, ids_control, ids_x;
    uint32_t int_sub;
    string word;

    // --- read the gate parameters ---
    read_reg_int(istr, ids_target);
    
    istr >> word;
    int_sub = get_value_from_word(word);

    read_reg_int(istr, ids_carry);

    // --- read the end of the gate structure ---
    YVVIv ids_control_it, ids_x_it;
    read_end_gate(istr, ids_control, ids_x, ids_control_it, ids_x_it);

    // --- gate ---
    x(ids_x); 
    comparator_fixed(ids_target, ids_carry, int_sub, ids_control, flag_inv);
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


void QCircuit::read_structure_gate_fourier(YISS istr, YCS path_in, YCB flag_inv)
{
    YVIv ids_target, ids_control, ids_x;

    // --- read target qubits ---
    read_reg_int(istr, ids_target);

    // --- read end of gate structure ---
    YVVIv ids_control_it, ids_x_it;
    read_end_gate(istr, ids_control, ids_x, ids_control_it, ids_x_it);

    // add the quantum Fourier circuit:
    x(ids_x);
    quantum_fourier(ids_target, ids_control, flag_inv, true);
    x(ids_x);
}


void QCircuit::read_structure_sin(YISS istr, YCS path_in, YCB flag_inv)
{
    YVIv ids_a, ids_main, ids_control, ids_x;
    qreal alpha_0, alpha;
    YVVIv ids_control_it, ids_x_it;
    string word;

    // --- read an ancilla qubit to put rotations there ---
    read_reg_int(istr, ids_a);

    // --- read angles ---
    istr >> word;
    alpha_0 = get_value_from_word(word);

    istr >> word;
    alpha = get_value_from_word(word);

    // --- read the condition qubits ---
    read_reg_int(istr, ids_main);

    // --- read end of gate structure ---
    read_end_gate(istr, ids_control, ids_x, ids_control_it, ids_x_it);

    // add the quantum Fourier circuit:
    x(ids_x);
    gate_sin(ids_a, ids_main, alpha_0, alpha, ids_control, flag_inv);
    x(ids_x);
}


void QCircuit::read_structure_gate_phase_estimation(YISS istr, YCS path_in, std::map<std::string, YSQ>& ocs, YCB flag_inv)
{
    YVIv ids_ta; // target qubits of the operator A, whose eigenphase we seek for;
    YVIv ids_ty;
    YVIv ids_cs, ids_x; 
    string name_A, name_INIT;

    // --- read target qubits of A ---
    read_reg_int(istr, ids_ta);

    // --- Read names of the operator A and INIT to find them among available circuits "ocs" ---
    try
    {
        istr >> name_A >> name_INIT;
    }
    catch(YCS e)
    {
        throw "PE definition: wrong format of the circuit names of the operators A and INIT: " + e;
    }

    // --- Read target qubit where the phase is to be written to ---
    read_reg_int(istr, ids_ty);
    
    // --- read the end of the gate structure ---
    YVVIv ids_control_it, ids_x_it;
    read_end_gate(istr, ids_cs, ids_x, ids_control_it, ids_x_it);

    // --- Find the appropriate circuits ---
    if(ocs.find(name_A) == ocs.end())
        throw string("PE definition: a circuit with the name ["s + name_A + "] is not found."s);
    if(ocs.find(name_INIT) == ocs.end())
        throw string("PE definition: a circuit with the name ["s + name_INIT + "] is not found."s);
    YSQ oc_A = ocs[name_A];
    YSQ oc_INIT = ocs[name_INIT];

    // --- add the phase estimation circuit ---
    x(ids_x);
    phase_estimation(ids_ta, oc_A, oc_INIT, ids_ty, ids_cs, flag_inv);
    x(ids_x); 
}


void QCircuit::read_structure_gate_qsvt(
    YISS istr, 
    YCS path_in, 
    std::map<std::string, YSQ>& ocs, 
    YCB flag_inv, 
    std::map<std::string, QSVT_pars>& map_qsvt_data
){
    string name_circuit;
    string be_name;
    YVIv ids_a_qsvt, ids_be;
    YVVIv ids_control_it, ids_x_it;
    YVIv ids_cs, ids_x;

    QSVT_pars data;

    // --- read QSVT parameters ---
    istr >> name_circuit;
    if(YMIX::compare_strings(name_circuit, unique_gates_names_))
        throw string("QSVT gate: the name ["s + name_circuit + "] has alredy been used for another gate."s);

    qsvt_read_parameters(name_circuit, data);

    // --- read QSVT ancilla ---
    read_reg_int(istr, ids_a_qsvt);
    if(YMIX::compare_strings(
        data.type, YVSv{"matrix-inversion", "gaussian-arcsin"}
    ))
        if(ids_a_qsvt.size() != 1)
            throw "QSVT circuit of type ["s + data.type + "] must have only a single ancilla specific qubit."s;
    if(YMIX::compare_strings(data.type, "hamiltonian-sim"))
        if(ids_a_qsvt.size() != 2)
            throw "QSVT circuit of type ["s + data.type + "] must have two ancilla specific qubits."s; 

    // --- read the name of the block-encoding oracle ---
    istr >> be_name;
    if(ocs.find(be_name) == ocs.end())
        throw string("QSVT definition: a circuit with the name ["s + be_name + "] is not found."s);
    YSQ oc_be = make_shared<QCircuit>(ocs[be_name]);

    // --- read qubits where the BE oracle must be placed ---
    read_reg_int(istr, ids_be);
    if(oc_be->get_n_qubits() != ids_be.size())
    {
        stringstream temp_sstr;
        temp_sstr << "QSVT definition: Number of qubits in the oracle is " << oc_be->get_n_qubits() << ", ";
        temp_sstr << "but it is indicated that the oracle is placed on " << ids_be.size() << " qubits.\n";
        throw string(temp_sstr.str());
    }

    // --- read the end of the gate structure ---
    read_end_gate(istr, ids_cs, ids_x, ids_control_it, ids_x_it);

    // --- add the QSVT circuit ---
    x(ids_x);
    if(data.parity == 0)
        qsvt_def_parity(data.angles_phis_even, ids_a_qsvt[0], ids_be, oc_be, ids_cs, flag_inv);
    if(data.parity == 1)
        qsvt_def_parity(data.angles_phis_odd, ids_a_qsvt[0], ids_be, oc_be, ids_cs, flag_inv);
    x(ids_x); 

    // store the QSVT data:
    map_qsvt_data[name_circuit] = data;
}


YQCP QCircuit::adder_1(YCVI ts, YCVI cs, YCB flag_inv)
{
    YVIv ids_target = YVIv(ts);
    uint32_t nt;

    if(!flag_inv)
    {
        // put the high-priority qubits at the beginning
        sort(ids_target.begin(), ids_target.end(), greater<int>());
        nt = ids_target.size();

        // add CNOT and X gates with control nodes
        for(unsigned i = 0; i < nt-1; ++i)
        {
            YVIv ids_cnot_cs = YVIv(ids_target.begin() + i + 1, ids_target.end());
            ids_cnot_cs.insert(ids_cnot_cs.end(), cs.begin(), cs.end());
            x(ids_target[i], ids_cnot_cs);
        }
        x(ids_target.back(), cs);
    }
    else
    {
        subtractor_1(ts, cs, false);
    }
    return get_the_circuit();
}


YQCP QCircuit::adder_2(YCVI ts, YCVI cs, YCB flag_inv)
{
    YVIv ids_target_init = YVIv(ts);
    uint32_t nt;

    if(!flag_inv)
    {
        // put the high-priority qubits at the beginning
        sort(ids_target_init.begin(), ids_target_init.end(), greater<int>());
        nt = ids_target_init.size()-1;
        YVIv ids_target(nt);
        copy(ids_target_init.begin(), ids_target_init.end()-1, ids_target.begin());
        
        // add CNOT and X gates with control nodes
        for(unsigned i = 0; i < nt-1; ++i)
        {
            YVIv ids_cnot_cs = YVIv(ids_target.begin() + i + 1, ids_target.end());
            ids_cnot_cs.insert(ids_cnot_cs.end(), cs.begin(), cs.end());
            x(ids_target[i], ids_cnot_cs);
        }
        x(ids_target.back(), cs);
    }
    else
    {
        subtractor_2(ts, cs, false);
    }
    return get_the_circuit();
}


YQCP QCircuit::adder_3(YCVI ts, YCVI cs, YCB flag_inv)
{
    YVIv ids_target = YVIv(ts);
    int32_t nt;
    int32_t sh = 2;

    if(!flag_inv)
    {
        // put the high-priority qubits at the beginning
        sort(ids_target.begin(), ids_target.end(), greater<int>());
        nt = ids_target.size();

        // add CNOT and X gates with control nodes
        for(int i = 0; i < (nt-1-sh); ++i)
        {
            YVIv ids_cnot_cs = YVIv(ids_target.begin() + i + 1, ids_target.end()-sh);
            ids_cnot_cs.insert(ids_cnot_cs.end(), cs.begin(), cs.end());
            x(ids_target[i], ids_cnot_cs);
        }
        if((nt-1-sh) >= 0) x(ids_target[nt-1-sh], cs);

        x(ids_target.back(), cs);
        for(int i = nt-2; i >= 0; --i)
        {
            YVIv ids_cnot_cs = YVIv(ids_target.begin() + i + 1, ids_target.end());
            ids_cnot_cs.insert(ids_cnot_cs.end(), cs.begin(), cs.end());
            x(ids_target[i], ids_cnot_cs);
        }
    }
    else
    {
        subtractor_3(ts, cs, false);
    }
    return get_the_circuit();
}


YQCP QCircuit::subtractor_1(YCVI ts, YCVI cs, YCB flag_inv)
{
    YVIv ids_target = YVIv(ts);
    uint32_t nt;

    if(!flag_inv)
    {
        // put the low-priority qubits at the beginning 
        sort(ids_target.begin(), ids_target.end());
        nt = ids_target.size();

        // add CNOT and X gates with control nodes 
        x(ids_target[0], cs);
        for(unsigned i = 1; i < nt; ++i)
        {
            YVIv ids_cnot_cs = YVIv(ids_target.begin(), ids_target.begin() + i);
            ids_cnot_cs.insert(ids_cnot_cs.end(), cs.begin(), cs.end());
            x(ids_target[i], ids_cnot_cs);
        }
    }
    else
    {
        adder_1(ts, cs, false);
    }
    return get_the_circuit();
}


YQCP QCircuit::subtractor_2(YCVI ts, YCVI cs, YCB flag_inv)
{
    YVIv ids_target_init = YVIv(ts);
    uint32_t nt;

    if(!flag_inv)
    {
        // put the low-priority qubits at the beginning 
        sort(ids_target_init.begin(), ids_target_init.end());
        nt = ids_target_init.size() - 1;
        YVIv ids_target(nt);
        copy(ids_target_init.begin()+1, ids_target_init.end(), ids_target.begin());

        // add CNOT and X gates with control nodes 
        x(ids_target[0], cs);
        for(unsigned i = 1; i < nt; ++i)
        {
            YVIv ids_cnot_cs = YVIv(ids_target.begin(), ids_target.begin() + i);
            ids_cnot_cs.insert(ids_cnot_cs.end(), cs.begin(), cs.end());
            x(ids_target[i], ids_cnot_cs);
        }
    }
    else
    {
        adder_2(ts, cs, false);
    }
    return get_the_circuit();
}


YQCP QCircuit::subtractor_3(YCVI ts, YCVI cs, YCB flag_inv)
{
    YVIv ids_target = YVIv(ts);
    uint32_t nt;
    uint32_t sh = 2;

    if(!flag_inv)
    {
        // put the low-priority qubits at the beginning 
        sort(ids_target.begin(), ids_target.end());
        nt = ids_target.size();

        for(unsigned i = nt-1; i > 0; --i)
        {
            YVIv ids_cnot_cs = YVIv(ids_target.begin(), ids_target.begin() + i);
            ids_cnot_cs.insert(ids_cnot_cs.end(), cs.begin(), cs.end());
            x(ids_target[i], ids_cnot_cs);
        }
        x(ids_target[0], cs);

        if(nt > sh) x(ids_target[sh], cs);
        for(unsigned i = sh+1; i < nt; ++i)
        {
            YVIv ids_cnot_cs = YVIv(ids_target.begin() + sh, ids_target.begin() + i);
            ids_cnot_cs.insert(ids_cnot_cs.end(), cs.begin(), cs.end());
            x(ids_target[i], ids_cnot_cs);
        }
    }
    else
    {
        adder_3(ts, cs, false);
    }
    return get_the_circuit();
}


YQCP QCircuit::adder(YCVI ts1, YCVI ts2, YCVI ts3, YCVI cs, YCB flag_inv, YCB flag_box)
{
    uint32_t nt = ts1.size();
    string oracle_name_tex = "ADD";
    YVIv ts_total; 
    ts_total.insert(ts_total.end(), ts1.begin(), ts1.end());
    ts_total.insert(ts_total.end(), ts2.begin(), ts2.end());
    ts_total.insert(ts_total.end(), ts3.begin()+1, ts3.end());
    ts_total.insert(ts_total.end(), ts3.begin(),   ts3.begin()+1);
    uint32_t nt_total = ts_total.size();

    // --- create an envelop circuit for the adder operator ---
    auto oc_adder = make_shared<QCircuit>(
        "ADD", env_, path_to_output_, nt_total
    );
    auto rc = oc_adder->add_register("rc", nt);
    auto rb = oc_adder->add_register("rb", nt);
    auto ra = oc_adder->add_register("ra", nt);

    for(uint32_t ii = 0; ii < nt; ii++)
    {
        oc_adder->x(rc[ii], YVIv {rb[ii], ra[ii]});
        oc_adder->x(rb[ii], YVIv {ra[ii]});
        if(ii > 0)
            oc_adder->x(rc[ii], YVIv {rc[ii-1], rb[ii]});
    } 

    if(nt > 1) 
    {
        oc_adder->x(rb[nt-1], YVIv {rc[nt-2]});
        for(uint32_t ii = nt-2; ii > 0; ii--)
        {
            oc_adder->x(rc[ii], YVIv {rc[ii-1], rb[ii]});
            oc_adder->x(rb[ii], YVIv {ra[ii]});
            oc_adder->x(rc[ii], YVIv {rb[ii], ra[ii]});
            oc_adder->x(rb[ii], YVIv {ra[ii]});
            oc_adder->x(rb[ii], YVIv {rc[ii-1]});
        }
        oc_adder->x(rb[0], YVIv {ra[0]});
        oc_adder->x(rc[0], YVIv {rb[0], ra[0]});
        oc_adder->x(rb[0], YVIv {ra[0]});
    }

    // --- invert the circuit if necessary ---
    if(flag_inv)
    {
        oc_adder->conjugate_transpose();
        oracle_name_tex += "^\\dagger";
    }

    // --- copy the adder circuit to the current circuit ---
    auto box = YSB(nullptr);
    if(flag_box)
        box = YMBo("ADD", ts_total, YVIv{}, oracle_name_tex);
    copy_gates_from(
        oc_adder,
        ts_total,
        box, 
        false,    
        cs        // add the control nodes at the cs qubits;
    );
    return get_the_circuit();
}


YQCP QCircuit::subtractor(YCVI ts1, YCVI ts2, YCVI ts3, YCVI cs, YCB flag_inv)
{
    auto cs_total_for_x = YVIv(cs);
    cs_total_for_x.push_back(ts3[0]);

    // if(flag_inv) 
    //     for(int ii = 0; ii < ts2.size(); ii++)
    //         x(ts2[ii], cs_total_for_x);

    x(ts2, cs);
    adder(ts1, ts2, ts3, cs, flag_inv);
    x(ts2, cs);

    // if(!flag_inv) 
    //     for(int ii = 0; ii < ts2.size(); ii++)
    //         x(ts2[ii], cs_total_for_x);

    return get_the_circuit();
}


YQCP QCircuit::adder_qft(YCVI qs_v1, YCVI qs_v2, YCI id_carry, YCVI cs, YCB flag_inv, YCB flag_box)
{
    int nv       = qs_v1.size();
    int nt_total = nv + 1;
    int n_circ_total = nt_total + nv;
    string oracle_name_tex = "ADDQFT";

    // --- create an envelop circuit for the adder operator ---
    auto oc_adder = make_shared<QCircuit>(
        "ADD-QFT", env_, path_to_output_, n_circ_total
    );
    auto loc_v2    = oc_adder->add_register("a", nv);
    auto loc_carry = oc_adder->add_register("c", 1)[0];
    auto loc_v1    = oc_adder->add_register("r", nv);
    auto loc_res = YVIv(loc_v1);
    loc_res.push_back(loc_carry);

    oc_adder->quantum_fourier(loc_res, YVIv{}, false, true);

    // starting from the most significant qubit (which is the carry bit):
    int n_phase_gates_per_qubit = 0;
    for(int i_qubit = nt_total-1; i_qubit >= 0; i_qubit--)
    {
        n_phase_gates_per_qubit++;
        for(int j_ph_gate = 0; j_ph_gate < n_phase_gates_per_qubit; j_ph_gate++)
        {
            if(j_ph_gate < nv)
            {
                double phase_curr = 2*M_PI / (1<<(n_phase_gates_per_qubit - j_ph_gate));
                oc_adder->phase(i_qubit, phase_curr, YVIv{loc_v2[j_ph_gate]});
            }
        }
    }
    oc_adder->quantum_fourier(loc_res, YVIv{}, true, true);

    // --- invert the circuit if necessary ---
    if(flag_inv)
    {
        oc_adder->conjugate_transpose();
        oracle_name_tex += "^\\dagger";
    }

    // --- copy the adder circuit to the current circuit ---
    auto ts_total = YVIv(qs_v1);
    ts_total.push_back(id_carry);
    ts_total.insert(ts_total.end(), qs_v2.begin(), qs_v2.end());

    auto box = YSB(nullptr);
    if(flag_box)
        box = YMBo("ADDFIXED", ts_total, YVIv{}, oracle_name_tex);
    copy_gates_from(
        oc_adder,
        ts_total,
        box, 
        false,    
        cs 
    );
    return get_the_circuit();
}


YQCP QCircuit::subtractor_qft(YCVI qs_v1, YCVI qs_v2, YCI id_sign, YCVI cs, YCB flag_inv, YCB flag_box)
{
    auto cs_total_for_x = YVIv(cs);
    cs_total_for_x.push_back(id_sign);

    // if(flag_inv) 
    //     for(int ii = 0; ii < qs_v1.size(); ii++)
    //         x(qs_v1[ii], cs_total_for_x);

    x(qs_v1, cs);
    adder_qft(qs_v1, qs_v2, id_sign, cs, flag_inv, flag_box);
    x(qs_v1, cs);

    // if(!flag_inv) 
    //     for(int ii = 0; ii < qs_v1.size(); ii++)
    //         x(qs_v1[ii], cs_total_for_x);

    return get_the_circuit();
}


YQCP QCircuit::adder_fixed(
        YCVI ids_target, YCI id_carry, YCU int_sub, 
        YCVI cs, YCB flag_inv, YCB flag_box
){
    int nv       = ids_target.size();
    int nt_total = nv + 1;
    string oracle_name_tex = "ADDFIXED";

    // --- represent the unisgned integer as a bitstring ---
    YVshv bitstring_int(nt_total);
    YMATH::intToBinary(int_sub, bitstring_int);

    // cout << "bitstring: ";
    // for(int i = 0; i < bitstring_int.size(); i++)
    //     cout << bitstring_int[i] << " ";
    // cout << endl;

    // --- create an envelop circuit for the adder operator ---
    auto oc_adder = make_shared<QCircuit>(
        "ADD-FIXED", env_, path_to_output_, nt_total
    );
    auto r_all = oc_adder->add_register("r", nt_total);

    oc_adder->quantum_fourier(r_all);

    // starting from the most significant qubit (which is the carry bit):
    int n_phase_gates_per_qubit = 0;
    for(int i_qubit = nt_total-1; i_qubit >= 0; i_qubit--)
    {
        n_phase_gates_per_qubit++;
        for(int j_ph_gate = 0; j_ph_gate < n_phase_gates_per_qubit; j_ph_gate++)
            if(bitstring_int[nt_total - j_ph_gate - 1] == 1)
            {
                double phase_curr = 2*M_PI / (1<<(n_phase_gates_per_qubit - j_ph_gate));
                oc_adder->phase(i_qubit, phase_curr);
            }
    }
    oc_adder->quantum_fourier(r_all, YVIv{}, true);

    // --- invert the circuit if necessary ---
    if(flag_inv)
    {
        oc_adder->conjugate_transpose();
        oracle_name_tex += "^\\dagger";
    }

    // --- copy the adder circuit to the current circuit ---
    auto ts_total = YVIv(ids_target);
    ts_total.push_back(id_carry);

    auto box = YSB(nullptr);
    if(flag_box)
        box = YMBo("ADDFIXED", ts_total, YVIv{}, oracle_name_tex);
    copy_gates_from(
        oc_adder,
        ts_total,
        box, 
        false,    
        cs 
    );
    return get_the_circuit();
}


YQCP QCircuit::subtractor_fixed(
        YCVI ids_target, YCI id_carry, YCU int_sub, 
        YCVI cs, YCB flag_inv
){
    auto ts_total = YVIv(ids_target);
    ts_total.push_back(id_carry);

    auto cs_total_for_x = YVIv(cs);
    cs_total_for_x.push_back(id_carry);

    // if(flag_inv) 
    //     for(int ii = 0; ii < ids_target.size(); ii++)
    //         x(ids_target[ii], cs_total_for_x);
            
    x(ts_total, cs);
    adder_fixed(ids_target, id_carry, int_sub, cs, flag_inv);
    x(ts_total, cs);

    // if(!flag_inv) 
    //     for(int ii = 0; ii < ids_target.size(); ii++)
    //         x(ids_target[ii], cs_total_for_x);

    return get_the_circuit();
}


YQCP QCircuit::comparator_fixed(
        YCVI ids_target, YCVI ids_carry, YCU int_sub, 
        YCVI cs, YCB flag_inv
){
    bool flag_inv_inv = !flag_inv;

    if(flag_inv_inv)
        cout << "comparator: here" << endl;

    auto cs_total_for_x = YVIv(cs);
    cs_total_for_x.push_back(ids_carry[0]);

    subtractor_fixed(ids_target, ids_carry[0], int_sub, cs);
    x(ids_carry[1], cs_total_for_x);
    subtractor_fixed(ids_target, ids_carry[0], int_sub, cs, true);
    return get_the_circuit();
}


YQCP QCircuit::swap(YCI t1, YCI t2, YVIv cs)
{ 
    YVIv ids_cs_1 = {t1};
    YVIv ids_cs_2 = {t2};
    ids_cs_1.insert(ids_cs_1.end(), cs.begin(), cs.end());
    ids_cs_2.insert(ids_cs_2.end(), cs.begin(), cs.end());
    return x(t2, ids_cs_1)->x(t1, ids_cs_2)->x(t2, ids_cs_1); 
}


YQCP QCircuit::quantum_fourier(YCVI ts, YCVI cs, YCB flag_inv, YCB flag_box)
{
    uint32_t nt = ts.size();
    uint32_t q1;
    int cq;
    qreal aa;
    string fourier_name_tex = "F";

    // --- create an envelop circuit for the Fourier operator ---
    auto oc_fourier = make_shared<QCircuit>(
        "F", env_, path_to_output_, nt
    );
    auto q = oc_fourier->add_register("q", nt);
    for(uint32_t ii = 0; ii < nt; ii++)
    {
        q1 = nt - 1 - ii;
        oc_fourier->h(q1);
        for(uint32_t jj = 1; jj < (nt - ii); jj++)
        {
            aa = M_PI / (1 << jj);
            oc_fourier->phase(q1, aa, YVIv{int(q1 - jj)});
        }
    }
    for(uint32_t ii = 0; ii < uint32_t(nt/2); ii++)
        oc_fourier->swap(q[ii], q[nt - 1 - ii]);

    // --- invert the circuit if necessary ---
    if(flag_inv)
    {
        oc_fourier->conjugate_transpose();
        fourier_name_tex += "^\\dagger";
    }

    // --- copy Fourier env. circuit to the current circuit ---
    auto box = YSB(nullptr);
    if(flag_box)
        box = YMBo("F", ts, YVIv{}, fourier_name_tex);
    copy_gates_from(
        oc_fourier,
        ts,
        box, 
        false,    
        cs        // add the control on the cs qubits;
    );
    return get_the_circuit();
}


YQCP QCircuit::gate_sin(
        YCVI anc, 
        YCVI conds, 
        YCQR alpha_0, 
        YCQR alpha, 
        YCVI cs, 
        YCB flag_inv, 
        YCB flag_box
){
    auto na = 1;
    auto n_cond = conds.size();
    auto n_tot = n_cond + na;
    string name_tex = "SIN";

    // --- all target qubits ---
    vector<int> qubits_tot = YVIv(conds);
    qubits_tot.insert(qubits_tot.end(), anc.begin(), anc.end());

    // --- create an envelop circuit for the sin gate ---
    auto oc_sin = make_shared<QCircuit>(
        "SIN", env_, path_to_output_, n_tot
    );
    auto a_loc    = oc_sin->add_register("a", na)[0];
    auto cond_loc = oc_sin->add_register("cond", n_cond);

    // cout << "\n qubits_tot:\n";
    // for(int ii = 0; ii < qubits_tot.size(); ii++)
    //     cout << qubits_tot[ii] << " ";

    // cout << "\n cond_loc:\n";
    // for(int ii = 0; ii < cond_loc.size(); ii++)
    //     cout << cond_loc[ii] << " ";

    // cout << "\n a_loc: " << a_loc << "\n";

    oc_sin->ry(a_loc, 2*alpha_0);
    for(auto ii = 0; ii < n_cond; ii++)
    {
        qreal aa = 2*alpha / pow(2., n_cond - 1 - ii);
        oc_sin->ry(a_loc, aa, YVIv{cond_loc[ii]});
    }
    oc_sin->x(a_loc);

    // --- invert the circuit if necessary ---
    if(flag_inv)
    {
        oc_sin->conjugate_transpose();
        name_tex += "^\\dagger";
    }

    // --- copy the env. circuit to the current circuit ---
    auto box = YSB(nullptr);
    if(flag_box)
        box = YMBo("SIN", qubits_tot, YVIv{}, name_tex);
    copy_gates_from(
        oc_sin,
        qubits_tot,
        box, 
        false,    
        cs        
    );
    return get_the_circuit();
}


YQCP QCircuit::phase_estimation(
    YCVI ta, 
    const std::shared_ptr<const QCircuit>& A, 
    const std::shared_ptr<const QCircuit>& INIT,
    YCVI ty, 
    YCVI cs, 
    YCB flag_inv,
    YCB flag_box
){
    string pe_name_tex = "PE";
    auto ny      = ty.size();
    auto nq_A    = A->get_n_qubits();
    if(ta.size() != nq_A)
    {
        string err_line;
        err_line  = "--- Error: setting the structure of the PE gate ---\n";
        err_line += "The subcircuit " + A->get_name() + " has " + to_string(nq_A) + " qubits, while" + 
            " one has indicated " + to_string(ta.size()) + " qubits to connect to.";
        throw err_line;
    }
    if(ta.size() != INIT->get_n_qubits())
    {
        string err_line;
        err_line  = "--- Error: setting the structure of the PE gate ---\n";
        err_line += "The subcircuit " + INIT->get_name() + " has " + to_string(INIT->get_n_qubits()) + " qubits, while" + 
            " one has indicated " + to_string(ta.size()) + " qubits to connect to.";
        throw err_line;
    }

    // --- create sequences from several A operators ---
    list<shared_ptr<QCircuit>> sequs_A;
    for(int iy = 0; iy < ny; iy++)
    {
        uint32_t N_rot = 1 << iy;
        auto c_temp = make_shared<QCircuit>(
            "A2^{"s+to_string(iy)+"}"s, env_, path_to_output_, nq_A
        );
        auto qa = c_temp->add_register("a", nq_A);
        for(int i_rot = 0; i_rot < N_rot; i_rot++)
            c_temp->copy_gates_from(A, qa);
        sequs_A.push_back(c_temp);
    }

    // --- create an envelop circuit for the phase estimation operator ---
    int count_c = -1;
    auto oc_pe = make_shared<QCircuit>(
        "PE", env_, path_to_output_, nq_A + ny
    );
    auto qy = oc_pe->add_register("y", ny);
    auto qa = oc_pe->add_register("a", nq_A);
    oc_pe->copy_gates_from(
        INIT, qa, 
        YMBo("INIT", qa, YVIv{})
        // YSB(nullptr)
    );
    oc_pe->h(qy);
    for(auto& one_sequ: sequs_A)
    {
        count_c++;
        auto ids_t_sequ = YVIv(qa);
        ids_t_sequ.push_back(qy[count_c]);
        oc_pe->copy_gates_from(
            one_sequ, 
            ids_t_sequ,  
            YMBo(one_sequ->get_name(), qa, YVIv{}, one_sequ->get_name()),
            // YSB(nullptr),
            false,
            {qy[count_c]}
        );
    }
    oc_pe->quantum_fourier(qy, YVIv{}, true, true);

    // --- invert the circuit if necessary ---
    if(flag_inv)
    {
        oc_pe->conjugate_transpose();
        pe_name_tex += "^\\dagger";
    }

    // --- copy the PE env. circuit to the current circuit ---
    auto circ_qubits = YVIv({ta});
    circ_qubits.insert(circ_qubits.end(), ty.begin(), ty.end());

    auto box = YSB(nullptr);
    if(flag_box)
        box = YMBo("PE", circ_qubits, YVIv{}, pe_name_tex);

    copy_gates_from(
        oc_pe,
        circ_qubits,
        box, 
        false,
        cs        // add the control on the cs qubits;
    );
    return get_the_circuit();
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
        throw "The constant with the name "s + const_name + " is not found."s;
    return constants_[const_name];
}


void  QCircuit::qsvt_read_parameters(YCS gate_name, QSVT_pars& data)
{
    string filename = gate_name + FORMAT_QSP;
    ifstream ff_qsvt(filename);
    if(!ff_qsvt.is_open()) throw "Error: there is not the file: "s + filename;

    string line, key_name;
    while (getline(ff_qsvt, line))
    {
        key_name = "";
        line = YMIX::remove_comment(line);
        if(line.find_first_not_of(' ') == string::npos)
            continue;

        istringstream iss(line);
        iss >> key_name;

        if(YMIX::compare_strings(key_name, "filename_angles"))
        {
            iss >> data.filename_angles;
            continue;
        }
    }
    ff_qsvt.close();

    if(data.filename_angles.empty())
    {
        throw "Name of the file with QSVT angles is not indicated in "s + filename;
        exit(-1);
    }
        
    // --- read the .hdf5 file with angles ---
    string temp = data.filename_angles;
    if(!YMIX::compare_strings(
        temp.substr(
            temp.size()-string(FORMAR_HDF5).size(),string(FORMAR_HDF5).size()
        ), FORMAR_HDF5)
    ) temp += FORMAR_HDF5;
    data.filename_angles = temp;

    YMIX::print_log("\nRead angles from the file: "s + data.filename_angles);

    YMIX::H5File ff;
    ff.set_name(data.filename_angles);
    ff.open_r();
    ff.read_scalar(data.type, "polynomial_type", "basic");
    ff.read_scalar(data.eps_qsvt, "eps", "basic");
    ff.read_scalar(data.f_par,    "par", "basic");

    if(YMIX::compare_strings(data.type, "matrix-inversion"))
    {
        ff.read_vector(data.angles_phis_odd, "odd", "angles");
        data.parity = 1;
    }
    if(YMIX::compare_strings(data.type, "gaussian-arcsin"))
    {
        ff.read_vector(data.angles_phis_even, "even", "angles");
        data.parity = 0;
    }
    if(YMIX::compare_strings(data.type, "hamiltonian-sim"))
    {
        ff.read_scalar(data.nt,                 "nt", "basic");
        ff.read_vector(data.angles_phis_odd,   "odd", "angles");
        ff.read_vector(data.angles_phis_even, "even", "angles");
        data.parity = -1;
    }
    ff.close();

    // --- print resulting parameters ---
    stringstream istr;
    istr << "--- QSVT gate with the name: " << gate_name << " ---\n";
    istr << "   QSVT type:  " << data.type     << ";\n";
    istr << "   QSVT error: " << data.eps_qsvt << ";\n";
    istr << "   Polynomial parity: "  << data.parity << ";\n";
    if(YMIX::compare_strings(data.type, "matrix-inversion"))
    {
        istr << "   number of angles: " << data.angles_phis_odd.size() << ";\n";
        istr << "   kappa: " << data.f_par << ";\n";
    }
    if(YMIX::compare_strings(data.type, "gaussian-arcsin"))
    {
        istr << "   number of angles: " << data.angles_phis_even.size() << ";\n";
        istr << "   mu: " << data.f_par << ";\n";
    }
    if(YMIX::compare_strings(data.type, "hamiltonian-sim"))
    {
        istr << "   number of angles: " << 
            data.angles_phis_odd.size() + data.angles_phis_even.size() << ";\n";
        istr << "   single time interval: " << data.f_par << ";\n";
        istr << "   number of the time intervals: " << data.nt << ";\n";
    }
    YMIX::print_log(istr.str());
}


YQCP QCircuit::qsvt_def_parity(
    YCVQ phis_in,
    YCI a_qsvt,
    YCVI qs_be_in, 
    const std::shared_ptr<const QCircuit> BE,
    YCVI cs, 
    YCB flag_inv,
    YCB flag_box
){
    YMIX::YTimer timer;
    auto phis = YVQv(phis_in);
    auto N_angles = phis.size();
    auto n_be     = BE->get_n_qubits();
    auto n_be_anc = BE->get_na();
    auto cs_total = YVIv(cs);

    string qsvt_name_tex = "QSVT";
    string be_box_name = "BE";
    string be_box_name_tex;
    string be_box_cc_name_tex;

    // N_angles = 4; /// for testing;

    // --- separate BE ancillae and input qubits ---
    auto qs_be     = YVIv(qs_be_in);
    vector<int> be_input(qs_be.begin(),                   qs_be.begin() + n_be - n_be_anc);
    vector<int>   be_anc(qs_be.begin() + n_be - n_be_anc, qs_be.end()                    );

    // --- pre-initialize the BE oracle ---
    // of the same size as the whole current circuit:
    auto oc_be = make_shared<QCircuit>(be_box_name, env_, path_to_output_, nq_);
    oc_be->add_register("r", nq_);
    oc_be->copy_gates_from(BE, qs_be, YSB(nullptr), flag_inv, cs);

    // --- create the complex-conjugated block-encoding oracle ---
    auto oc_be_inv = make_shared<QCircuit>(oc_be);
    oc_be_inv->conjugate_transpose();

    // --- QSVT circuit ---
    be_box_name_tex    = be_box_name;
    be_box_cc_name_tex = be_box_name + "^\\dagger"s;
    if(flag_inv)
    {
        reverse(phis.begin(), phis.end());
        be_box_name_tex    = be_box_cc_name_tex;
        be_box_cc_name_tex = be_box_name;
    }
    cs_total.insert(cs_total.end(), be_anc.begin(), be_anc.end());
    sort(cs_total.begin(), cs_total.end());

    timer.StartPrint("Creating the QSVT circuit... ");
    h(a_qsvt, cs);
    x(be_anc, cs);
    x(a_qsvt, cs_total)->rz(a_qsvt, 2*phis[0], cs, flag_inv)->x(a_qsvt, cs_total);
    x(be_anc, cs);
    for(uint32_t count_angle = 1; count_angle < N_angles; ++count_angle)
    {
        insert_gates_from(
            oc_be.get(), 
            make_shared<Box__>(be_box_name, qs_be, cs, be_box_name_tex)
            // YSB(nullptr)
        );

        if((N_angles % 2) == 0)
        {
            if(!flag_inv)
            {
                if(count_angle == (N_angles - 1)) 
                    z(a_qsvt, cs);
            }
            else
            {
                if(count_angle == 1) 
                    z(a_qsvt, cs);
            }
        }
        x(be_anc, cs);
        x(a_qsvt, cs_total)->rz(a_qsvt, 2*phis[count_angle], cs, flag_inv)->x(a_qsvt, cs_total);
        x(be_anc, cs);
        
        count_angle += 1;
        if(count_angle < N_angles)
        {
            insert_gates_from(
                oc_be_inv.get(), 
                make_shared<Box__>(be_box_name, qs_be, cs, be_box_cc_name_tex)
                // YSB(nullptr)
            );
            x(be_anc, cs);
            x(a_qsvt, cs_total)->rz(a_qsvt, 2*phis[count_angle], cs, flag_inv)->x(a_qsvt, cs_total);
            x(be_anc, cs);
        }
    }
    h(a_qsvt, cs);
    timer.StopPrint();
    
    return get_the_circuit();
}
















// YQCP QCircuit::qsvt_def_parity(
//     YCVQ phis,
//     YCVI a_qsvt,
//     YCVI qs_be, 
//     const std::shared_ptr<const QCircuit> BE,
//     YCVI cs, 
//     YCB flag_inv,
//     YCB flag_box
// ){
//     string qsvt_name_tex = "QSVT";
//     auto n_be = BE->get_n_qubits();
//     auto n_be_anc = BE->get_na();

//     // --- pre-initialize the BE oracle ---
//     auto oc_be = make_shared<QCircuit>(
//         "BE", env_, path_to_output_, n_be
//     );
//     auto be_qubits  = oc_be->add_register("qs-be", n_be);
//     oc_be->copy_gates_from(
//         BE, 
//         be_qubits,
//         make_shared<Box__>("U", be_qubits, YVIv {})
//     );

//     // --- create the complex-conjugated block-encoding oracle ---
//     auto oc_be_inv = make_shared<QCircuit>(oc_be);
//     oc_be_inv->conjugate_transpose();

//     // --- initialize the envelop circuit for the QSVT procedure ---
//     auto oc_qsvt = make_shared<QCircuit>(
//         "QSVT", env_, path_to_output_, n_be + 1
//     );
//     auto loc_a_qsvt = oc_qsvt->add_register("a-qsvt", 1);
//     auto loc_a_be   = oc_qsvt->add_register("a-be", n_be_anc); // ancillae of the BE oracle;
//     auto loc_in     = oc_qsvt->add_register("in", n_be - n_be_anc); // input (main) qubits of the BE oracle;

//     auto loc_all_be = YVIv(loc_in);
//     loc_all_be.insert(loc_all_be.end(), loc_a_be.begin(), loc_a_be.end());

//     int q = loc_a_qsvt[0];
    
//     // cout << "\nloc_a_be:\n";
//     // for(int ii = 0; ii < loc_a_be.size(); ii++)
//     //     cout << loc_a_be[ii] << " ";
    
//     // cout << "\nloc_in:\n";
//     // for(int ii = 0; ii < loc_in.size(); ii++)
//     //     cout << loc_in[ii] << " ";

//     // cout << "\nbe_qubits:\n";
//     // for(int ii = 0; ii < be_qubits.size(); ii++)
//     //     cout << be_qubits[ii] << " ";

//     // cout << "\n q: " << q << "\n";

//     // cout << "a_qsvt:\n";
//     // for(int ii = 0; ii < a_qsvt.size(); ii++)
//     //     cout << a_qsvt[ii] << " ";

//     // cout << "\nqs_be:\n";
//     // for(int ii = 0; ii < qs_be.size(); ii++)
//     //     cout << qs_be[ii] << " ";
//     // cout << "\n";

//     // --- form the QSVT circuit for the odd polynomial (even number of angles) ---
//     YMIX::YTimer timer;
//     qreal aa;
//     auto N_angles = phis.size();

//     // N_angles = 4;

//     timer.StartPrint("Creating the QSVT circuit... ");
//     oc_qsvt->h(q);

//     aa = 2*phis[0];
//     oc_qsvt->x(loc_a_be);
//     oc_qsvt->x(q, loc_a_be)->rz(q, aa)->x(q, loc_a_be);
//     oc_qsvt->x(loc_a_be);
//     for(uint32_t count_angle = 1; count_angle < N_angles; ++count_angle)
//     {
//         if(flag_inv)
//             oc_qsvt->copy_gates_from(
//                 oc_be, 
//                 loc_all_be,
//                 make_shared<Box__>("U", be_qubits, YVIv {})
//             );
//         else
//             oc_qsvt->insert_gates_from(
//                 oc_be.get(), 
//                 make_shared<Box__>("U", be_qubits, YVIv {})
//             );
        
//         if(count_angle == (N_angles - 1)) // set the Z-gate only if N_angles is even (odd polynomial)
//             oc_qsvt->z(q);
//         aa = 2*phis[count_angle];
//         oc_qsvt->x(loc_a_be);
//         oc_qsvt->x(q, loc_a_be)->rz(q, aa)->x(q, loc_a_be);
//         oc_qsvt->x(loc_a_be);
        
//         count_angle += 1;
//         if(count_angle < N_angles)
//         {
//             if(flag_inv)
//                 oc_qsvt->copy_gates_from(
//                     oc_be_inv,
//                     loc_all_be, 
//                     make_shared<Box__>("iU", be_qubits, YVIv {})
//                 );
//             else
//                 oc_qsvt->insert_gates_from(
//                     oc_be_inv.get(), 
//                     make_shared<Box__>("iU", be_qubits, YVIv {})
//                 );
//             aa = 2*phis[count_angle];
//             oc_qsvt->x(loc_a_be);
//             oc_qsvt->x(q, loc_a_be)->rz(q, aa)->x(q, loc_a_be);
//             oc_qsvt->x(loc_a_be);
//         }
//     }
//     oc_qsvt->h(q);
//     timer.StopPrint();
    
//     // --- invert the QSVT circuit if necessary ---
//     if(flag_inv)
//     {
//         timer.StartPrint("Inversion of the QSVT circuit... ");
//         YMIX::print_log("inversion of the QSVT circuit");
//         oc_qsvt->conjugate_transpose();
//         qsvt_name_tex += "^\\dagger";
//         timer.StopPrint();
//     }

//     // --- copy the QSVT env. circuit to the current circuit ---
//     timer.StartPrint("Transferring the QSVT circuit to the the circuit ["s + name_ + "]... ");
//     auto circ_qubits = YVIv(qs_be);
//     circ_qubits.insert(circ_qubits.end(), a_qsvt.begin(), a_qsvt.end());

//     auto box = YSB(nullptr);
//     if(flag_box)
//         box = YMBo("QSVT", circ_qubits, YVIv{}, qsvt_name_tex);
//         if(flag_inv)
//             box = YMBo("QSVTcc", circ_qubits, YVIv{}, qsvt_name_tex);

//     copy_gates_from(
//         oc_qsvt,
//         circ_qubits,
//         box, 
//         false,
//         cs  
//     );
//     timer.StopPrint();

//     return get_the_circuit();
// }