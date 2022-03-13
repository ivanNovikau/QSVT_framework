#ifndef ORACLETOOL_H
#define ORACLETOOL_H

#include "QCircuit.h"

/*
The idea of the OracleTool is to read a circuit from a file, launch it several times
to output results for different input states and see states at different circuit points. 
*/
class OracleTool__
{
public:
    /**
     * @brief Construct a new OracleTool__ object
     * @param project_name  is a project name that defines names of input and output files.
     * @param path_to_inputs is path to input files.
     * @param flag_compute_output if True, compute output states, otherwise only input states
     * @param flag_circuit print or not .circuit files;
     * @param flag_tex print or not .tex files;
     */
    OracleTool__(
        const QuESTEnv& env, 
        YCS project_name, 
        YCS path_to_inputs, 
        YCB flag_compute_output,
        YCB flag_circuit,
        YCB flag_tex
    );
    ~OracleTool__();

    void read_input_file(YS data);
    void read_circuit_structure_from_file(YCS data);

    /**
     * @brief Launch the circuit with different input states.
     */
    void launch();

    /**
     * @brief Return shared pointer to the oracle.
     */
    inline YSQ get_oracle_copy() const { return std::make_shared<QCircuit>(oc_to_launch_); }

private:
    void read_constants(YISS str);
    void read_circuit_declaration(YISS istr);
    void read_circuit_structure(YISS istr);
    void read_input_states(YISS istr);

    void read_gate(YISS istr, YPQC oc, YCB flag_inv=false);
    void read_subcircuit(YISS istr, YPQC oc, YCB flag_inv=false);

    void read_state(YISS istr);

    qreal get_value_from_word(YCS word);

private:
    YGlobalVariables gv; // global variables;
    QuESTEnv env_; 
    std::string pname_;       // project name
    std::string path_inputs_; // path to input files
    std::string ifname_;      // name of an input file

    unsigned nc_; // number of circuits
    std::map<std::string, YSQ> ocs_; // circuits
    std::map<std::string, qreal> constants_; // constants
    YSQ oc_to_launch_ = nullptr; // a circuit to launch

    bool flag_compute_output_; // if True, compute output states from an oracle;
    bool flag_circuit_; // print or not .circuit files;
    bool flag_tex_; // print or not .circuit files;

    /** 
     * several initial states:
     * every init. state is represented by several registers, where several qubits might be set to 1.
     */
    std::vector<std::map<std::string, std::vector<int>>> init_states_;


};


#endif

