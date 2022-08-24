#include "../include/circuittool.h"
using namespace std;


CircuitTool__::CircuitTool__(
    const QuESTEnv& env, 
    YCS pname, 
    YCS path_to_inputs, 
    YCB flag_random,
    YCB flag_compute_output,
    YCB flag_print_output,
    YCB flag_tex,
    YCB flag_hdf5
) : BaseTool__(
    env, pname, path_to_inputs, 
    flag_compute_output, flag_print_output, 
    flag_random, flag_tex, false, flag_hdf5
){
    format_file_ = FORMAT_CIRCUIT;
    flag_random_ = flag_random;

    if(!flag_random)
        read_data();
    else
        create_random_circuit();
}


CircuitTool__::~CircuitTool__()
{
    YMIX::print_log( "*** Destruction of the circuit tool is done. ***");
}


 void CircuitTool__::read_circuit_structure_from_file(YCS data)
 {
    istringstream istr(data);
    string word, line;
    string circuit_name;
    uint32_t nq;
    uint64_t counter_line = 0;

    // ----------------------------------------------------------------------------
    // --- Read the circuit name and the number of qubits ---
    counter_line++;
    {
        getline(istr, line);
        istringstream iss(line);
        iss >> circuit_name >> nq;
        oc_to_launch_ = make_shared<QCircuit>(
                circuit_name, env_, path_inputs_, nq, 
                map<string, qreal>(),
                flag_circuit_, flag_tex_, flag_layers_
            );
        YMIX::print_log(
            "Circuit [" + circuit_name + "] with " + to_string(nq) + " qubits is created."
        );
    }
    
    // ----------------------------------------------------------------------------
    // --- Read the description of the circuit registers ---
    counter_line++;
    {
        string reg_name;
        string line_log;
        uint32_t nq_in_reg;
        bool flag_anc;

        getline(istr, line);
        istringstream iss(line);
        iss >> word;
        if(!YMIX::compare_strings(word, "QubitRegisterNames"))
            throw string(
                "Information about the circuit registers is not found "s + 
                "(the keyword QubitRegisterNames is absent)."s
            );
        while(iss >> reg_name >> nq_in_reg >> flag_anc)
        {
            line_log = "add ";
            if(flag_anc) line_log+= "an ancilla "s;
            else line_log += "a "s;
            line_log += "register [" + reg_name + "] with " + to_string(nq_in_reg) + " qubits;";
            YMIX::print_log( line_log);
            oc_to_launch_->add_register(reg_name, nq_in_reg, flag_anc);
        }
        oc_to_launch_->save_regs();
        oc_to_launch_->set_standart_output_format();
    }

    // ----------------------------------------------------------------------------
    // --- Read gates ---
    string gate_name, box_name, line_conj, line_start_end;
    bool flag_conj, flag_start;
    int64_t id_layer, n;
    shared_ptr<Gate__> gate;
    while (getline(istr, line))
    {
        counter_line++;
        flag_conj = false;
        flag_start = false;

        // cout << "--- Line " << counter_line << " ---\n";

        istringstream iss(line);
        iss >> gate_name;

        if(YMIX::compare_strings(gate_name, "Box"))
        {
            iss >> line_start_end >> box_name; 
            if(YMIX::compare_strings(line_start_end, "start"))
                flag_start = true;
        }

        iss >> line_conj;
        if(YMIX::compare_strings(line_conj, "conj"))
            flag_conj = true;

        iss >> word;
        if(!YMIX::compare_strings(word, "layer"))
            throw string("Line " + to_string(counter_line) + ": the keyword layer is absent.");
        iss >> id_layer;

        iss >> word;
        if(!YMIX::compare_strings(word, "targets"))
            throw string("Line " + to_string(counter_line) + ": the keyword targets is absent.");
        iss >> n;
        vector<int> targets(n);
        for(auto i = 0; i < n; i++)
            iss >> targets[i];

        iss >> word;
        if(!YMIX::compare_strings(word, "controls"))
            throw string("Line " + to_string(counter_line) + ": the keyword controls is absent.");
        iss >> n;
        vector<int> controls(n);
        for(auto i = 0; i < n; i++)
            iss >> controls[i];

        iss >> word;
        if(!YMIX::compare_strings(word, "pars"))
            throw string("Line " + to_string(counter_line) + ": the keyword pars is absent.");
        iss >> n;
        vector<qreal> pars(n);
        for(auto i = 0; i < n; i++)
            iss >> pars[i];
            
        if(YMIX::compare_strings(gate_name, "Box"))
            gate = make_shared<Box__>(box_name, targets, YVIv {}, "");
        if(YMIX::compare_strings(gate_name, "X"))
            gate = make_shared<X__>(targets[0]);
        if(YMIX::compare_strings(gate_name, "Y"))
            gate = make_shared<Y__>(targets[0]);
        if(YMIX::compare_strings(gate_name, "Z"))
            gate = make_shared<Z__>(targets[0]);
        if(YMIX::compare_strings(gate_name, "H"))
            gate = make_shared<H__>(targets[0]);
        if(YMIX::compare_strings(gate_name, "Rx"))
            gate = make_shared<Rx__>(targets[0], pars[0]);
        if(YMIX::compare_strings(gate_name, "Ry"))
            gate = make_shared<Ry__>(targets[0], pars[0]);
        if(YMIX::compare_strings(gate_name, "Rz"))
            gate = make_shared<Rz__>(targets[0], pars[0]);
        if(YMIX::compare_strings(gate_name, "Phase"))
            gate = make_shared<Phase__>(targets[0], pars[0]);
        if(YMIX::compare_strings(gate_name, "Rc"))
            gate = make_shared<Rc__>(targets[0], pars[0], pars[1]);

        gate->add_control_qubits(controls);
        gate->set_flag_start(flag_start);
        gate->set_layer(id_layer);
        if(flag_conj) gate->conjugate_transpose();
        oc_to_launch_->add_gate(gate);
    }
 }


void CircuitTool__::create_random_circuit()
{
    uint32_t nq = 0;
    uint64_t n_max_gates = 0;
    string circuit_name = pname_;

    // --- READ .random file to now the desired shape of the random circuit ---
    string file_name = path_inputs_ + "/" + pname_ + FORMAT_RANDOM;
    string data;
    read_input_file(data, file_name);

    istringstream istr(data);
    istr >> nq >> n_max_gates;

    // --- CONSTRUCT the random circuit ---
    {
        oc_to_launch_ = make_shared<QCircuit>(
                    circuit_name, env_, path_inputs_, nq, 
                    map<string, qreal>(),
                    flag_circuit_, flag_tex_, flag_layers_
                );
        YMIX::print_log(
                "Circuit [" + circuit_name + "] with " + to_string(nq) + " qubits is created."
            );
        oc_to_launch_->add_register("r", nq, false);
        oc_to_launch_->save_regs();
        oc_to_launch_->set_standart_output_format();
    }

    uint32_t n_avail_gates = Gate__::avail_gate_names_.size();
    srand (time(NULL));
    uint32_t tq1, tq2, cq;
    string gate_name;
    qreal rot_angle;
    for(uint64_t counter_gate = 0; counter_gate < n_max_gates; counter_gate++)
    {
        int temp = rand() % n_avail_gates;
        gate_name = Gate__::avail_gate_names_[temp];
        tq1 = rand() % nq;
        rot_angle = float(rand())/float((RAND_MAX)) * (2*M_PI);

        if(YMIX::compare_strings(gate_name, X__::name_shared_))
            oc_to_launch_->x(tq1);
        if(YMIX::compare_strings(gate_name, Y__::name_shared_))
            oc_to_launch_->y(tq1);
        if(YMIX::compare_strings(gate_name, Z__::name_shared_))
            oc_to_launch_->z(tq1);
        if(YMIX::compare_strings(gate_name, H__::name_shared_))
            oc_to_launch_->h(tq1);
        if(YMIX::compare_strings(gate_name, Rx__::name_shared_))
            oc_to_launch_->rx(tq1, rot_angle);
        if(YMIX::compare_strings(gate_name, Ry__::name_shared_))
            oc_to_launch_->ry(tq1, rot_angle);
        if(YMIX::compare_strings(gate_name, Rz__::name_shared_))
            oc_to_launch_->rz(tq1, rot_angle);
        if(YMIX::compare_strings(gate_name, Phase__::name_shared_))
            oc_to_launch_->phase(tq1, rot_angle);
        if(YMIX::compare_strings(gate_name, "SWAP"))
        {
            tq2 = rand() % nq;
            while(tq2 == tq1) tq2 = rand() % nq;
            oc_to_launch_->swap(tq1, tq2);
        }
        if(YMIX::compare_strings(gate_name, "CNOT"))
        {
            cq = rand() % nq;
            while(cq == tq1) cq = rand() % nq;
            oc_to_launch_->x(tq1, YVIv {int(cq)});
        }
        if(YMIX::compare_strings(gate_name, "CH"))
        {
            cq = rand() % nq;
            while(cq == tq1) cq = rand() % nq;
            oc_to_launch_->h(tq1, YVIv {int(cq)});
        }
        if(YMIX::compare_strings(gate_name, "CRy"))
        {
            cq = rand() % nq;
            while(cq == tq1) cq = rand() % nq;
            oc_to_launch_->ry(tq1, rot_angle, YVIv {int(cq)});
        }
    }
}


void CircuitTool__::launch()
{
    // print the gates:
    oc_to_launch_->print_gates();

    // ---- store basic data:
    if(flag_hdf5_)
    {
        hfo_.open_w();

        // number of qubits
        hfo_.add_scalar(oc_to_launch_->get_n_qubits(), "nq", "basic");
        hfo_.add_scalar(oc_to_launch_->get_na(), "na", "basic");

        // register names
        string res_lin = "";
        for(auto const& reg_name: oc_to_launch_->get_reg_names())
            res_lin += reg_name + ", ";
        res_lin.pop_back(); res_lin.pop_back();
        hfo_.add_scalar(res_lin, "register-names", "basic");

        // number of qubits in every register:
        hfo_.add_vector(oc_to_launch_->get_standart_output_format(), "register-nq", "basic");

        // number of initial states:
        hfo_.add_scalar(1, "n-init-states", "states");

        hfo_.close();
    }

    // calculate the circuit states:
    if(flag_compute_output_)
    {
        YMIX::YTimer timer;
        YMIX::StateVectorOut outF;

        // --- For consistency, store the initial state, which is always  |000..00> ---
        if(flag_hdf5_)
        {
            oc_to_launch_->get_state(outF);
            hfo_.open_w();
            hfo_.add_vector(outF.ampls,  "initial-amplitudes-0"s, "states");
            hfo_.add_matrix(outF.states, "initial-states-0"s,     "states");
            hfo_.close(); 
        }

        // --- Calculate the output state ---
        YMIX::print_log("\n--- Calculating the output state... ---");
        timer.Start();
        oc_to_launch_->generate();
        oc_to_launch_->get_state(outF);
        timer.Stop();
        YMIX::print_log("Done: duration: " + timer.get_dur_str_s());
        if(flag_print_output_) YMIX::print_log("Output state:\n" + outF.str_wv);

        // --- Store the output state ---
        if(flag_hdf5_)
        {
            hfo_.open_w();
            hfo_.add_vector(outF.ampls,  "output-all-amplitudes-0"s, "states");
            hfo_.add_matrix(outF.states, "output-all-states-0"s,     "states");
            hfo_.close(); 
        }
    }
}