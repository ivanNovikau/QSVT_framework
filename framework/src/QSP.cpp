#include "../include/QSP.h"

using namespace std;


QSP__::QSP__(const QuESTEnv& env, YCS project_name, YCS path_inputs)
    : QSVT__(env, project_name, path_inputs)
{
    nq_ = 2; // the circuit has two specific ancillar (+ qubits from the oracle)
    type_ = "QSP";
}


QSP__::~QSP__()
{
    // destroy inner circuits:
    cw_.reset();
    icw_.reset();
    cu_.reset();
    icu_.reset();
}


bool QSP__::read_special_parameters(std::istringstream& iss, std::string& key_name)
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
            YMIX::print_log(env_, "WARNING: nt_points is non-positive, set to 1;");
        }
        if(dt < 0)
        {
            dt = abs(dt);
            YMIX::print_log(env_, "WARNING: dt is negative, set to |dt|;");
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
        cf << "\n--- QSP simulation of the time evolution ---\n";
        cf << "dt = " << f_par_ << "\n";
        cf << "Time grid: [" << t_grid_[0] << ", " << t_grid_.back() << "] with " << t_grid_.size() 
            << " points;\n";
        cf << "-------------------------------------------------------\n";

        return true;
    }
    return false;
}


void QSP__::read_angles()
{
    YMIX::print_log(env_, "Reading QSP angles... ");
    if(env_.rank == 0)
    {
        string fname_full;
        string fname_res;
        ifstream ff;
        qreal f_par, eps;

        // calculate angles for Hamiltonian simulation:
        fname_res = 
            "angles_t" + to_string(int(1000*f_par_)) + 
            "_eps" + to_string(int(-log10(eps_))) + 
            FORMAT_ANGLES;
        fname_full = path_inputs_ + "/" + fname_res;

        YMIX::print_log(env_, "Angles are read from the file "s + fname_full);
        ff.open(fname_full.c_str());
        if(ff.is_open())
        {
            ff >> f_par >> eps >> N_angles_;
            angles_phis_ = YVQv(N_angles_);
            for(unsigned i = 0; i < N_angles_; ++i)
                ff >> angles_phis_[i];
            ff.close();
        }
        else
        {
            YMIX::print_log_err("There is not the following file\n" + fname_full);
        }

        // --- recheck the parameters ---
        if(!YMATH::is_zero(abs(eps - eps_)))
        {
            throw string("eps does not correspond to the QSVD angles.");
        }
        if(!YMATH::is_zero(abs(f_par - f_par_)))
        {
            throw string("QSVD parameter (time step or kappa or etc.) does not correspond to the QSVD angles.");
        }

        ostringstream ocs;
        ocs << "function parameter = " << f_par << ";   ";
        ocs << "eps = "            << eps << ";   ";
        ocs << "N of angles = "    << N_angles_ << ";\n";
        YMIX::print_log(env_, "\n"s + ocs.str());
    }

    // --- Save the angles to .hdf5 file ---
    YMIX::print_log(env_, "\nWriting angles to .hdf5 file");
    hfo_.open_w();
    hfo_.add_vector(angles_phis_, "qsp-angles-one-time-step", "basic");
    hfo_.close();
}


void QSP__::save_basic_specific_data()
{
    hfo_.add_scalar(f_par_, "dt", "basic");
    hfo_.add_vector(t_grid_, "time-grid", "basic");
    hfo_.add_scalar(eps_, "qsp-initial-precision", "basic");
}


void QSP__::create_circuit_components()
{
    // --- Create the controlled oracle ---
    auto cu = make_shared<QCircuit>(u_);
    cu->controlled(regs_["qb"][0]);
    cu_ = make_shared<const QCircuit>(cu);

    // create the inverse controlled oracle:
    auto conjCU = make_shared<QCircuit>(cu_);
    conjCU->conjugate_transpose();
    icu_ = make_shared<QCircuit>(conjCU);

    YMIX::print_log(env_, "Qubitisation...");
    create_W();
    create_iW();
}


void QSP__::create_W()  // original W
{ 
    auto a_ancs = u_->get_ancillae();
    YVIv aq_ancs = YVIv(a_ancs);
    aq_ancs.push_back(regs_["qb"][0]);

    // --- create an iterate W ---
    YMIX::print_log(env_, "-> cW creation");
    cw_ = make_shared<QCircuit>(
        "W", env_, path_inputs_, nq_,
        map<string, qreal>(),
        flag_circuit_, false, flag_layers_
    );

    int  rb   = cw_->add_register("b", 1)[0];
    int  rq   = cw_->add_register("q", 1)[0];
    YVIv ureg = cw_->add_register("ureg", nq_-2); 
    cw_->save_regs();

    auto ts_box = YMATH::get_range(0,nq_-2);
    YVIv cs_box = {rq};
    
    // add cU 
    YMIX::print_log(env_, "Adding cU...");
    cw_->x(rq);
    cw_->copy_gates_from(cu_, YMATH::get_range(0, nq_), make_shared<Box__>("U", ts_box, cs_box));
    cw_->x(rq);

    // add cU*
    YMIX::print_log(env_, "Adding cU*...");
    auto conjU = make_shared<QCircuit>(cu_);
    conjU->conjugate_transpose();
    cw_->copy_gates_from(conjU, YMATH::get_range(0, nq_), make_shared<Box__>("U*", ts_box, cs_box));

    // add a reflector operator 
    cw_->h(rq);

    for(auto& anc1: aq_ancs)
        cw_->x(anc1);
    cw_->phase(rq, M_PI, a_ancs);
    for(auto& anc1: aq_ancs)
        cw_->x(anc1);

    cw_->x(rq);
    cw_->phase(rq, M_PI);
    cw_->x(rq);

    cw_->h(rq);

    // make the oracle W controlled 
    cw_->controlled({rb});

    cw_->print_gates();
}


void QSP__::create_iW()
{
    YMIX::print_log(env_, "-> icW creation...\n");
    icw_ = make_shared<QCircuit>(cw_, "W*");
    icw_->conjugate_transpose();
    icw_->print_gates();
}


void QSP__::simulation()
{
    YMIX::YTimer timer;
    unsigned nt = t_grid_.size();
    int b = nq_ - 1;
    int q = nq_ - 2;
    auto ts_box = YMATH::get_range(0,nq_-1);
    YVIv cs_box = {b};
    qreal aa;
    YMIX::StateVectorOut outZ;

    // iterate on time steps
    for(unsigned i = 0; i < nt; ++i)
    {
        if(env_.rank == 0)
        {
            ostringstream ocs;
            ocs << "\n--- Time step " << i << " ---\n";
            ocs << "   t = "          << t_grid_[i] << ";";
            YMIX::print_log(env_, ocs.str());
        }

        // structure a QSP circuit:
        timer.StartPrint(env_, "Structuring a QSP circuit... ");
        oc_->h(b)->h(q);
        for(unsigned count_angle = 0; count_angle < int((N_angles_-1)/2); ++count_angle)
        {
            aa = angles_phis_[2*count_angle];
            oc_->rz(b, -aa)->h(b)->x(b);
            oc_->insert_gates_from(
                cw_.get(), 
                make_shared<Box__>("W", ts_box, cs_box)
            );
            oc_->phase(b, -M_PI_2);
            oc_->x(b)->h(b)->rz(b, aa);

            // oc_->x(0); // for non-hermitian


            aa = angles_phis_[2*count_angle+1];
            oc_->rz(b, -aa)->h(b);
            oc_->insert_gates_from(
                icw_.get(),
                make_shared<Box__>("W*", ts_box, cs_box)
            );
            oc_->phase(b, M_PI_2);
            oc_->h(b)->rz(b, aa);

            // oc_->x(0); // for non-hermitian
        }
        aa = angles_phis_[N_angles_-1];
        oc_->rz(b, aa);
        oc_->h(b)->h(q);

        // oc_->x(0); // for non-hermitian

        timer.StopPrint(env_);

        // --- Generate the circuit ---
        timer.StartPrint(env_, "Calculating the QSP circuit... ");
        oc_->generate();
        oc_->get_state(outZ, true);
        timer.StopPrint(env_);

        if(flag_print_zero_states_)
        {
            cout << "resulting zero-ancilla state ->\n";
            cout << outZ.str_wv << endl;
        }

        // --- Save the resulting states to .hdf5 file ---
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