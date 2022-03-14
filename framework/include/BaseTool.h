#ifndef BASETOOL_H
#define BASETOOL_H

#include "QCircuit.h"


class BaseTool__
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
    BaseTool__(
        const QuESTEnv& env, 
        YCS project_name, 
        YCS path_to_inputs, 
        YCB flag_compute_output,
        YCB flag_circuit,
        YCB flag_tex,
        YCB flag_layer
    );
    ~BaseTool__();

    inline YSQ get_circuit() const { return std::make_shared<QCircuit>(oc_to_launch_); }
    virtual void launch() = 0;
    
protected:
    void read_data();
    void read_input_file(YS data);
    virtual void read_circuit_structure_from_file(YCS data) = 0;
    
protected:
    QuESTEnv env_; 
    std::string pname_;       // project name
    std::string path_inputs_; // path to input files
    std::string ifname_;      // name of the input file
    std::string format_file_; // format of the input file;

    unsigned nc_; // number of circuits
    std::map<std::string, YSQ> ocs_; // circuits
    YSQ oc_to_launch_ = nullptr; // a circuit to launch

    bool flag_compute_output_; // if True, compute output states from an oracle;
    bool flag_circuit_; // print or not .circuit files;
    bool flag_tex_; // print or not .circuit files;
    bool flag_layers_; // to calculate or not the layers;
};

#endif