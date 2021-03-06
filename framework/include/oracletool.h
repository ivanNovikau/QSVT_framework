#ifndef ORACLETOOL_H
#define ORACLETOOL_H

#include "BaseTool.h"

class OracleTool__ : public BaseTool__
{
public:
    /**
     * @brief Construct a new OracleTool__ object
     * @param[in] project_name  is a project name that defines names of input and output files.
     * @param[in] path_to_inputs is path to input files.
     * @param[in] flag_compute_output if True, compute output states, otherwise only input states
     * @param[in] flag_circuit print or not .circuit files;
     * @param[in] flag_tex print or not .tex files;
     * @param[in] flag_layers to calculate or not the layer for each gate;
     */
    OracleTool__(
        const QuESTEnv& env, 
        YCS project_name, 
        YCS path_to_inputs, 
        YCB flag_compute_output,
        YCB flag_circuit,
        YCB flag_tex,
        YCB flag_layer
    );
    ~OracleTool__();
    void launch();

protected:
    void read_circuit_structure_from_file(YCS data);
    

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
    // dictuinary of constants to create the oracle:
    std::map<std::string, qreal> constants_; 

    // several initial states:
    // every init. state is represented by several registers, 
    //  where several qubits might be set to 1.
    std::vector<std::map<std::string, std::vector<int>>> init_states_;
};
#endif

