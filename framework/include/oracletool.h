#ifndef ORACLETOOL_H
#define ORACLETOOL_H

#include "BaseTool.h"

class OracleTool__ : public BaseTool__
{
public:
    /**
     * @param[in] project_name  is a project name that defines names of input and output files.
     * @param[in] path_to_inputs is path to input files.
    //  * @param[in] flag_compute_output if True, compute output states, otherwise only input states
    //  * @param[in] flag_print_output if True, print output states on screen;
    //  * @param[in] flag_circuit print or not .circuit files;
    //  * @param[in] flag_tex print or not .tex files;
    //  * @param[in] flag_layers to calculate or not the layer for each gate;
    //  * @param[in] flag_hdf5  to create the hdf5 file;
     */
    OracleTool__(
        const QuESTEnv& env, 
        YCS project_name, 
        YCS path_to_inputs 
    );
    ~OracleTool__();
    void launch();

protected:
    void read_circuit_structure_from_file(YCS data);
    

private:
    void read_constants(YISS str);
    void read_options(YISS str);
    void read_circuit_declaration(YISS istr);
    void read_circuit_structure(YISS istr);
    void read_main_circuit(YISS istr);
    void read_gate(YISS istr, YPQC oc, YCB flag_inv=false);
    void read_subcircuit(YISS istr, YPQC oc, YCB flag_inv=false);
    void read_state(YISS istr);
    void read_state_init_file();
    qreal get_value_from_word(YCS word);
    void  calc(std::shared_ptr<QCircuit>& u_work, YCI count_init_state);

private:
    // dictionary with constants to create the oracle:
    std::map<std::string, qreal> constants_; 

    // several initial states:
    // every init. state is represented by several registers, 
    //  where several qubits might be set to 1.
    std::vector<std::map<std::string, std::vector<int>>> init_states_;

    // amplitudes of the initial state read from the initial state file:
    std::vector<qreal> init_ampl_vec_real_;
    std::vector<qreal> init_ampl_vec_imag_;

    // true if the initial state is read from the .init_state file:
    bool flag_init_state_file_;

    /**
     * none          - do not compute/print any output states (no printing on the sceen);
     * all           - compute all output states (and separately, zero-ancilla states);
     *                     however, do not print separately the zero-ancilla states;
     * zero-ancillae - compute/print states where all ancillae are in the zero state
     * (remark, ancilla qubits are assumed to have higher priority than the main qubits).
     */ 
    std::string sel_compute_output_; 
    std::string sel_print_output_; 

    /**
     * To calculate probabilities of all states for the given qubits, \p focus_qubits_.
     */
    bool flag_prob_;
    std::vector<int> focus_qubits_;

    // QSVT data
    std::map<std::string, QSVT_pars> qsvt_data_;


};
#endif

