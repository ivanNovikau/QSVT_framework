#include "../include/QDYN.h"

using namespace std;


QDYN__::QDYN__(const QuESTEnv& env, YCS project_name, YCS path_inputs)
    : QSVT__(env, project_name, path_inputs)
{
    nq_ = 2; // the circuit has two specific ancillar (+ qubits from the oracle)
    type_ = "QSVT";
}


bool QDYN__::read_special_parameters(std::istringstream& iss, std::string& key_name)
{
    // description of the time grid:
    if(YMIX::compare_strings(key_name, "t_grid"))
    {
        int nt_points;
        qreal dt;
        iss >> dt >> nt_points;
        
        if(nt_points <= 0)
        {
            nt_points = 1;
            YMIX::print_log( "WARNING: nt_points is non-positive, set to 1;");
        }
        if(dt < 0)
        {
            dt = abs(dt);
            YMIX::print_log( "WARNING: dt is negative, set to |dt|;");
        }
        if(YMATH::is_zero(dt))
        {
            YMIX::print_log_err("Error: dt == 0;");
        }

        f_par_ = dt;
        t_grid_ = YVQv(nt_points);
        for(unsigned i = 1; i <= nt_points; ++i)  
            t_grid_[i-1] = dt * i;

        YMIX::LogFile cf;
        cf << "\n--- QSVD simulation of the time evolution ---\n";
        cf << "dt = " << f_par_ << "\n";
        cf << "Time grid: [" << t_grid_[0] << ", " << t_grid_.back() << "] with " << t_grid_.size() 
            << " points;\n";
        cf << "-------------------------------------------------------\n";

        return true;
    }
    return false;
}


void QDYN__::read_angles()
{
    YMIX::print_log( "Reading QSVD angles... ");
    if(env_.rank == 0)
    {
        read_angles_def_parity("even", angles_phis_even_, N_angles_even_);
        read_angles_def_parity("odd",  angles_phis_odd_,  N_angles_odd_);
    }
}


void QDYN__::save_basic_specific_data()
{
    hfo_.add_scalar(f_par_, "dt", "basic");
    hfo_.add_vector(t_grid_, "time-grid", "basic");
}


void QDYN__::create_circuit_components()
{
    // cos
    create_circuit_component_def_parity("even", angles_phis_even_, N_angles_even_, ceven_, false);
    create_controlled_component(ceven_, ceven_controlled_);

    // i*sin
    create_circuit_component_def_parity("odd", angles_phis_odd_, N_angles_odd_, codd_, false);
    create_controlled_component(codd_, codd_controlled_);
}


void QDYN__::simulation()
{
    YMIX::YTimer timer;
    unsigned nt = t_grid_.size();
    int b = nq_ - 1;
    auto ts_box = YMATH::get_range(0,nq_-2);
    YVIv cs_box = {b};
    YMIX::StateVectorOut outZ;
    for(unsigned i = 0; i < nt; ++i)
    {
        if(env_.rank == 0)
        {
            ostringstream ocs;
            ocs << "\n--- Time step " << i << " ---\n";
            ocs << "   t = "          << t_grid_[i] << ";";
            YMIX::print_log( ocs.str());
        }

        // structure the QSVT circuit:
        timer.StartPrint("Structuring the QSVT circuit... ");

        // oc_->x(b);  // to have exp(-i*H*t)
        oc_->h(b);

        oc_->x(b);
        oc_->insert_gates_from(
            ceven_controlled_.get(), 
            make_shared<Box__>("E", ts_box, cs_box)
        );
        oc_->x(b);

        oc_->insert_gates_from(
            codd_controlled_.get(), 
            make_shared<Box__>("O", ts_box, cs_box)
        );

        oc_->h(b);
        timer.StopPrint();

        // --- generate the circuit ---
        timer.StartPrint("Calculating the QSVT circuit... ");
        oc_->generate();
        oc_->get_state(outZ, true);
        timer.StopPrint();

        if(flag_print_zero_states_)
        {
            cout << "resulting zero-ancilla state ->\n";
            cout << outZ.str_wv << endl;
        }

        // --- Save the resulting zero-ancilla states to .hdf5 file ---
        hfo_.open_w();
        hfo_.add_vector(outZ.ampls,  "t-step-" + to_string(i) + "--output-amplitudes", "states");
        hfo_.add_matrix(outZ.states, "t-step-" + to_string(i) + "--output-states", "states");
        hfo_.close(); 

        // --- Get the full state ---
        if(flag_print_all_states_)
        {
            YMIX::StateVectorOut outF;
            oc_->get_state(outF);
            cout << "resulting full state ->\n";
            cout << outF.str_wv << endl;
        }
    }
    
}