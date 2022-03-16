#include "../include/circuittool.h"
using namespace std;


CircuitTool__::CircuitTool__(
    const QuESTEnv& env, 
    YCS pname, 
    YCS path_to_inputs, 
    YCB flag_compute_output,
    YCB flag_tex
) : BaseTool__(
    env, pname, path_to_inputs, 
    flag_compute_output, false, flag_tex, false
){
    format_file_ = FORMAT_CIRCUIT;
    read_data();
}

CircuitTool__::~CircuitTool__()
{
    YMIX::print_log(env_, "*** Destruction of the circuit tool is done. ***");
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
            env_, 
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
            YMIX::print_log(env_, line_log);
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


 void CircuitTool__::launch()
 {
     // print the gates:
     oc_to_launch_->print_gates();

     // calculate the output state:
     if(flag_compute_output_)
     {
        YMIX::YTimer timer_gen;
        string str_wv, str_wv_nz;

        YMIX::print_log(env_, "\n--- Calculating the output state... ---");
        timer_gen.Start();
        oc_to_launch_->generate();
        oc_to_launch_->get_wavefunction(
            oc_to_launch_->get_standart_output_format(), str_wv, str_wv_nz, 3
        );
        timer_gen.Stop();
        YMIX::print_log(env_, "Done: duration: " + timer_gen.get_dur_str_s());
        YMIX::print_log(env_, "Output state:\n" + str_wv_nz);
     }
 }