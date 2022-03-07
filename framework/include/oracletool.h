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
     * @param flag_test is a flag, if True, output additional information from the program
     * @param flag_compute_iterator_output compute an output from an iterator
     */
    OracleTool__(
        const QuESTEnv& env, 
        YCS project_name, 
        YCS path_to_inputs, 
        const bool& flag_compute_output = true,
        const bool& flag_test = false,
        const bool& flag_compute_iterator_output = false
    );
    ~OracleTool__();

    void read_input_file(YS data);
    void read_circuit_structure_from_file(YCS data);

    /**
     * @brief Launch the circuit with different input states.
     */
    void launch();

    /**
     * @brief Return a shared pointer to an oracletool.
     */
    inline YSQ get_copy_of_an_oracle() const { return std::make_shared<QCircuit>(oc_to_launch_); }

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

    bool flag_test_; // flag to output more data for the code analysis
    bool flag_compute_output_; // if True, compute output states from an oracle
    bool flag_compute_iterator_output_; // if True, compute output states from an iterator

    /** several initial states:
     * every init. state represents several registers, where several qubits might be set to 1.
    */
    std::vector<std::map<std::string, std::vector<int>>> init_states_;


};


#endif

