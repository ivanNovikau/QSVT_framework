#include "../include/oracletool.h"
using namespace std;

/**
 * @brief To launch the oracletool.
 * @param argv oracletool [project_name] [path_to_input_files] [flag_compute_output] 
 *  [flag-circuit] [flag-tex] [tex-circuit-length] [flag-layers]
 */
int main(int argc, char *argv[])
{
    QuESTEnv env = createQuESTEnv();

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
        cout << "Project name: " << argv[1] << endl;
        cout << "Path to input files: " << argv[2] << "/" << endl;
    }

    // project name
    string pname(argv[1]);
    
    // path to input files:
    string path_input(argv[2]);
    if(env.rank == 0)
    {
        YMIX::LogFile::name_global_ = path_input + "/" + pname + ".clog";
        YMIX::LogFile cf(true);
    }
    YMIX::print_log(env, "Number of ranks: " + to_string(env.numRanks));

    // flag to compute output from the oracle
    bool flag_compute_output = true;
    if(argc > 3) flag_compute_output = stoi(string (argv[3]));
    if(flag_compute_output) YMIX::print_log(env, "Flag compute output oracle = true");
    else                    YMIX::print_log(env, "Flag compute output oracle = false");

    // flag print .circuit files:
    bool flag_circuit = true;
    if(argc > 4) flag_circuit = stoi(string (argv[4]));
    if(flag_circuit) YMIX::print_log(env, "Print .circuit files.");
    else             YMIX::print_log(env, "Do not print the .circuit files.");

    // flag print .tex files:
    bool flag_tex = true;
    if(argc > 5) flag_tex = stoi(string (argv[5]));
    if(flag_tex) YMIX::print_log(env, "Print .tex files.");
    else         YMIX::print_log(env, "Do not print the .tex files.");

    // length of the circuit in the .tex file:
    if(argc > 6) 
        YGV::tex_circuit_length = stoi(string (argv[6]));
    if(flag_tex) 
        YMIX::print_log(
            env, 
            "Circuit length in the .tex files = " + to_string(YGV::tex_circuit_length)
        );  

    // flag to calculate layers:
    bool flag_layers = false;
    if(argc > 7) flag_layers = stoi(string (argv[7]));
    if(flag_layers) YMIX::print_log(env, "Calculate layers for each gate.");
    else            YMIX::print_log(env, "Do not calculate layers.");
    
    try
    {
        // create the tool for the oracle analysis: 
        OracleTool__ oo = OracleTool__(
            env, 
            pname, 
            path_input, 
            flag_compute_output,
            flag_circuit,
            flag_tex,
            flag_layers,
            true
        );
        
        // launch the circuit with different input states:
        oo.launch();
    }
    catch(YCS e)
    {
        if(env.rank == 0) std::cerr << "\n" << e << endl;
        destroyQuESTEnv(env);
        return -1;
    }
    catch(const std::exception& e)
    {
        std::cerr << "General error:\n" << e.what() << '\n';
        destroyQuESTEnv(env);
        return -1;
    }
    
    // free the memory
    destroyQuESTEnv(env);
    
    return 0;
}