#include "../include/QMI.h"

using namespace std;


QMI__::QMI__(const QuESTEnv& env, YCS project_name, YCS path_inputs)
    : QSVT__(env, project_name, path_inputs)
{
    nq_ = 1; // the circuit has a single specific ancilla (+ qubits from the oracle)
    type_ = "QSVT";
}



bool QMI__::read_special_parameters(std::istringstream& iss, std::string& key_name)
{
    // read the condition number:
    if(YMIX::compare_strings(key_name, "kappa"))
    {    
        iss >> f_par_;

        if(f_par_ < 0)
        {
            YMIX::print_log_err(env_, "Error: kappa is negative;");
        }

        if(YMATH::is_zero(f_par_))
        {
            YMIX::print_log_err(env_, "Error: kappa == 0;");
        }

        YMIX::LogFile cf;
        cf << "\n--- Matrix inversion ---\n";
        cf << "condition number: " << f_par_ << ";\n";
        cf << "-------------------------------------------------------\n";

        return true;
    }
    return false;
}

void QMI__::read_angles()
{
    YMIX::print_log(env_, "Reading QSVD angles... ");
    if(env_.rank == 0)
    {
        read_angles_def_parity("odd", angles_phis_odd_, N_angles_odd_);
    }
}


void QMI__::save_basic_specific_data()
{
    hfo_.add_scalar(f_par_, "kappa", "basic");
}


void QMI__::create_circuit_components()
{
    create_circuit_component_def_parity("odd", angles_phis_odd_, N_angles_odd_, codd_, false);
}


void QMI__::simulation()
{
    YMIX::YTimer timer;
    auto ts_box = YMATH::get_range(0,nq_-1);
    YMIX::StateVectorOut outZ;

    // --- Resulting QSVT circuit ---
    timer.StartPrint(env_, "Structuring the QSVT circuit... ");
    oc_->insert_gates_from(
        codd_.get(), 
        make_shared<Box__>("QC", ts_box, YVIv {})
    );
    timer.StopPrint(env_);

    // --- generate the circuit ---
    timer.StartPrint(env_, "Calculating the QSVT circuit... ");
    oc_->generate();
    oc_->get_state(outZ, true);
    timer.StopPrint(env_);
    if(flag_print_zero_states_)
    {
        cout << "resulting zero-ancilla state ->\n";
        cout << outZ.str_wv << endl;
    }
    
    // --- Save the resulting zero-ancilla states to .hdf5 file ---
    hfo_.open_w();
    hfo_.add_vector(outZ.ampls, "output-amplitudes", "states");
    hfo_.add_matrix(outZ.states, "output-states", "states");
    hfo_.close(); 

    // --- Print full state ---
    if(flag_print_all_states_)
    {
        YMIX::StateVectorOut outF;
        oc_->get_state(outF);
        cout << "resulting full state ->\n";
        cout << outF.str_wv << endl;
    } 
}