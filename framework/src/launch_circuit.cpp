#include "../include/circuittool.h"
using namespace std;

/**
 * @brief To launch the qc_circuit program.
 * @param argv qc_circuit [project_name] [path_to_input_files] [flag-output] [flag-tex] [tex-circuit-length]
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
        cout << "Circuit name: " << argv[1] << endl;
        cout << "Path to input files: " << argv[2] << "/" << endl;
    }

    // --------------------------------------------------------
    uint32_t id_arg;

    // project name
    id_arg = 1;
    string pname(argv[id_arg]);
    
    // path to input files:
    id_arg += 1;
    string path_input(argv[id_arg]);
    if(env.rank == 0)
    {
        YMIX::LogFile::name_global_ = path_input + "/" + pname + ".clog";
        YMIX::LogFile cf(true);
    }
    YMIX::print_log(env, "Number of ranks: " + to_string(env.numRanks));
    YMIX::print_log(env, "\n--- Initial flags ---");

    // flag to compute output from the oracle
    id_arg += 1;
    bool flag_compute_output = true;
    if(argc > id_arg) flag_compute_output = stoi(string (argv[id_arg]));
    if(flag_compute_output) YMIX::print_log(env, "-> print output state.");
    else                    YMIX::print_log(env, "-> do not print output state.");

    // flag print .tex files:
    id_arg += 1;
    bool flag_tex = true;
    if(argc > id_arg) flag_tex = stoi(string (argv[id_arg]));
    if(flag_tex) YMIX::print_log(env, "-> print the .tex file.");
    else         YMIX::print_log(env, "-> do not print the .tex file.");

    // length of the circuit in the .tex file:
    id_arg += 1;
    if(argc > id_arg) 
        YGV::tex_circuit_length = stoi(string (argv[id_arg]));
    if(flag_tex) 
        YMIX::print_log(
            env, 
            "-> circuit length in the .tex file = " + to_string(YGV::tex_circuit_length)
        );  

    YMIX::print_log(env, "---\n");
    // --------------------------------------------------------
    try
    {
        CircuitTool__ oo = CircuitTool__(
            env, 
            pname, 
            path_input, 
            flag_compute_output,
            flag_tex
        );
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
    destroyQuESTEnv(env);
    return 0;
}