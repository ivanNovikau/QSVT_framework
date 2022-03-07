#include "../include/QDYN.h"
#include "../include/QMI.h"
#include "../include/QSP.h"
#include "../include/oracletool.h"
using namespace std;

/**
 * @brief To launch the QSVD framework.
 * @param argv qsp, [project_name], [path_to_input_files]
 */
int main(int argc, char *argv[])
{
    QuESTEnv env = createQuESTEnv();

    // read [project_name], [path_to_input_files]:
    if(env.rank == 0)
    {
        if(argc < 2)
        {
            cerr << "Project name is missing" << endl;
            exit(-1);
        }
        if(argc < 3)
        {
            cerr << "Path to input files is missing" << endl;
            exit(-1);
        }
        if(argc < 4)
        {
            cerr << "It is not indicated what to simulate: qsvt or qsp" << endl;
            exit(-1);
        }
        cout << "Project name: "        << argv[1] << endl;
        cout << "Path to input files: " << argv[2] << "/" << endl;
        cout << "Case to simulate: "    << argv[3] << endl;
    }

    // --- project name ---
    string pname(argv[1]);
    
    // --- path to input files ---
    string path_input(argv[2]);

    // --- case to simulate ---
    string case_to_sim(argv[3]); // qsp, qsvt-dyn, qsvt-mi

    // --- create a QSP log file ---
    if(env.rank == 0)
    {
        // create
        YMIX::LogFile::name_global_ = path_input + "/" + pname + FORMAT_LOG;
        YMIX::LogFile cf(true);
    }

    // --- write an initial information about the project ---
    string str_date_time;
    YMIX::get_current_date_time(str_date_time);
    YMIX::print_log(env, str_date_time);
    YMIX::print_log(env, "Number of ranks: " + to_string(env.numRanks));
    YMIX::print_log(env, "Case to simulate: " + case_to_sim);


    // --------------------------------------------------------------------------------
    // --- QSVD computation ---
    try
    {
        // --- Launch a tool to create an oracle from a given description ---
        // -> oracle description is given in a [path_input]/[pname].oracle file;
        OracleTool__* oo = new OracleTool__(
            env, 
            pname, 
            path_input, 
            false, 
            false
        );
        YSQ oc_oracle = oo->get_copy_of_an_oracle();
        delete oo;

        // if(case_to_sim.compare("qsp") == 0)
        // {
        //     cout << "*** QSP SIMULATION ***" << endl;


        //     // --- Create a QSP framework ---
        //     QSP__ qsp = QSP__(env, pname, path_input);

        //     // --- pass the oracle to the QSP framework to create a QSP circuit ---
        //     qsp.create_circuit("QSP", oc_oracle);

        //     // --- Print initial QSP parameters ---
        //     qsp.print_init_data();

        //     // --- Set an initial state ---
        //     qsp.set_init_vector();

        //     // --- Qubitisation ---
        //     qsp.qubitization();

        //     // // --- test oracle ---
        //     // qsp.test_U(oc_oracle.get());

        //     // // --- test iterator ---
        //     // qsp.test_W("W");

        //     // --- get measurement ---
        //     // set_measurement(qsp);

        //     // --- perform the QSP simulation ---
        //     qsp.QuantumSignalProcessing();
        // }

        QSVT__* qsvd;

        // --- QSP: Simulation of the time evolution using QSP ---
        if(case_to_sim.compare("qsp") == 0)
        {
            qsvd = new QSP__(env, pname, path_input);
        }

        // --- QSVT: Simulation of the time evolution ---
        if(case_to_sim.compare("qsvt-dyn") == 0)
        {
            qsvd = new QDYN__(env, pname, path_input);
        }
        
        // -- QSVT: Simulation of the matrix inversion ---
        if(case_to_sim.compare("qsvt-mi") == 0)
        {
            qsvd = new QMI__(env, pname, path_input);
        }

        // --- Read QSVD angles ---
        qsvd->read_angles();

        // --- pass the oracle to the QSVD framework to create a QSP circuit ---
        qsvd->create_circuit("QSVT", oc_oracle);

        // --- Print initial QSVT parameters ---
        qsvd->print_init_data();

        // --- Save basic data to the .hdf5 file ---
        qsvd->save_basic_data();

        // --- Set the initial state ---
        qsvd->set_init_vector();

        // --- perform the QSVD simulation ---
        qsvd->launch();

        delete qsvd;
        
    }
    catch(const string& e)
    {
        if(env.rank == 0) std::cerr << "\n" << "Error: " << e << endl;
        return -1;
    }
    catch(const std::exception& e)
    {
        std::cerr << "General error:\n" << e.what() << '\n';
        return -1;
    }
    // catch failure caused by the H5File operations
    catch( H5::FileIException error )
    {
       error.printError();
       return -1;
    }
    // catch failure caused by the DataSet operations
    catch( H5::DataSetIException error )
    {
       error.printError();
       return -1;
    }
    // catch failure caused by the DataSpace operations
    catch( H5::DataSpaceIException error )
    {
       error.printError();
       return -1;
    }
    // catch failure caused by the DataSpace operations
    catch( H5::DataTypeIException error )
    {
       error.printError();
       return -1;
    }

    // free the memory
    destroyQuESTEnv(env);

    return 0;
}