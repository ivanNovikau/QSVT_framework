#ifndef QSVT_H
#define QSVT_H

#include "../include/oracletool.h"

/**
 * @brief Quantum Singular Value Transformation (QSVT) framework.
 */
class QSVT__
{
public:
    /**
     * @brief Construct a Quantum Signal Processing object.
     * @param env QuEST environment.
     * @param fname_input name of an input file for with the QSP parameters.
     * @param path_inputs path to input files.
     */
    QSVT__(const QuESTEnv& env, YCS fname_input, YCS path_inputs);

    ~QSVT__();

    void read_main_parameters();
    void prepare_hdf5_files();
    void print_init_data();

    void read_block_encoding_oracle();

    void create_circuit();

    /**
     * @brief Read an initial state.
     */
    void set_init_vector();

    /**
     * @brief Insert gates from a circuit into the current circuit
     * without modification of their qubits and without copying them
     * @param[in] c circuit, from where gates are about to be appended to the current circuit.
     */
    void append(QCircuit* c);

    /**
     * @brief Calculate QSP angles at a current time point. 
     */
    

    /**
     * @brief Construct the final circuit for the quantum signal value decompos.
     * @param[in] state_to_output state to output;
     */
    void launch();

    /** 
     * @return Get a  number of ancilla qubits in the whle QSP circuit.
     */
    inline unsigned get_na() const { return oc_->get_na(); }

    inline std::string get_name_oracle() const { return u_->get_name(); }


    /**
     * @brief Save basic circuit parameters to the .hdf5 file.
     */
    void save_basic_data();

    /**
     * @brief Read QSVT aangles.
     */
    virtual void read_angles() = 0;


    void print_gates();

protected:
    void read_angles_def_parity(YCS line_parity, YVQ phis, uint32_t& N_angles);
    
    void create_circuit_component_def_parity(
        YCS line_parity, YCVQ phis, YCU N_angles, 
        YSQ& circ,  
        YCB flag_imaginary
    );
    void create_controlled_component(const YSQ circ, YSQ& circ_controlled);
    void save_restart_data();
    
    virtual void create_circuit_components() = 0;
    virtual bool read_special_parameters(std::istringstream& iss, std::string& key_name) = 0;
    virtual void save_basic_specific_data() = 0;
    virtual void simulation() = 0;

protected:
    QuESTEnv env_;
    YMIX::YTimer timer_;

    std::string type_;

    std::string project_name_;

    std::string path_inputs_; // path to input files
    std::string fname_input_; // input file name;
    std::string fname_init_;  // file with an initial state (represented as a vector);

    bool flag_restart_; // flag to read initial state from a restart file;

    bool flag_layers_; // to calculate or not the layers;

    // an .hdf5 file where output data are to write down:
    YMIX::H5File hfo_;

    // restart .hdf5 file: every next restart rewrites a restart file;
    // QSP reads previous restart file and at the end of simulation rewrites it:
    YMIX::H5File rf_;

    // indicates an index number of a restart:
    // 0 corresponds to the absence of previous restarts:
    unsigned restart_counter_; 

    qreal f_par_; // QMI: condition number of the matrix to inverse;
                // QDYN and QSP: length of one time step;
    qreal eps_; // error of the polynomial approximation;

    std::unique_ptr<QCircuit> oc_; // the QSP quantum circuit;
    int nq_; // total number of qubits in the circuit;

    std::map<std::string, YVIv> regs_; // registers;

    uint32_t N_angles_even_; // number of angles to build cos;
    YVQv angles_phis_even_; // angles to construct cos;

    uint32_t N_angles_odd_; // number of angles to build i*sin;
    YVQv angles_phis_odd_; // angles to construct i*sin;

    bool flag_circuit_; // print or not the QSVT (QSP) .circuit file;
    bool flag_print_zero_states_; // print or not zero states to the out log file;
    bool flag_print_all_states_; // print or not all states to the out log file;

    YSCQ oc_init_; // oracle to initialize the QSVT circuit;
    YSCQ u_; // oracle circuit;
    YSCQ iu_; // inversed oracle circuit;

    YSQ codd_; // odd component of the QSVT circuit;
    YSQ ceven_; // even component of the QSVT circuit;

    YSQ codd_controlled_; // controlled odd component of the QSVT circuit;
    YSQ ceven_controlled_; // controlled even component of the QSVT circuit;

    qreal coef_time_norm_; // time normalization factor;

    SEL_INIT_STATE_PREP sel_init_; // how to create the initial state;


};




#endif