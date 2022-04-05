#include "../include/circuittool.h"
using namespace std;

/**
 * @brief To launch the qc_circuit program.
 * @param argv qc_circuit [project_name] [path_to_input_files] 
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

    // --- project name ---
    id_arg = 1;
    string pname(argv[id_arg]);
    
    // --- path to input files ---
    id_arg += 1;
    string path_input(argv[id_arg]);
    if(env.rank == 0)
    {
        YMIX::LogFile::name_global_ = path_input + "/" + pname + ".clog";
        YMIX::LogFile cf(true);
    }
    YMIX::print_log(env, "Number of ranks: " + to_string(env.numRanks));
    
    // --------------------------------------------------------
    // --- Read flags ---
    bool flag_random = false;
    bool flag_compute_output = true;
    bool flag_print_output = true;
    bool flag_tex = true;

    id_arg += 1;
    while(id_arg < (argc - 1))
    {
        if(YMIX::compare_strings(argv[id_arg], "-flag_random"))
        {   
            id_arg += 1;
            flag_random = stoi(string (argv[id_arg]));
        }
        if(YMIX::compare_strings(argv[id_arg], "-flag_compute_output"))
        {
            id_arg += 1;
            flag_compute_output = stoi(string (argv[id_arg]));
        }
        if(YMIX::compare_strings(argv[id_arg], "-flag_print_output"))
        {
            id_arg += 1;
            flag_print_output = stoi(string (argv[id_arg]));
        }
        if(YMIX::compare_strings(argv[id_arg], "-flag_tex"))
        {
            id_arg += 1;
            flag_tex = stoi(string (argv[id_arg]));
        }
        if(YMIX::compare_strings(argv[id_arg], "-tex_CL"))
        {
            id_arg += 1;
            YGV::tex_circuit_length = stoi(string (argv[id_arg]));
        }

        id_arg += 1;
    }
    YMIX::print_log(env, "\n--- Initial flags ---");
    if(flag_random)          YMIX::print_log(env, "-> create a random circuit.");
    if(!flag_compute_output) YMIX::print_log(env, "-> do not compute output states.");
    if(!flag_print_output)   YMIX::print_log(env, "-> do not print output states on screen.");
    if(!flag_tex)            YMIX::print_log(env, "-> do not write the .tex file.");
    if(flag_tex) 
        YMIX::print_log(
            env, 
            "-> length of a single row in the .tex file = " + to_string(YGV::tex_circuit_length)
        );
    YMIX::print_log(env, "---\n");

    // --------------------------------------------------------
    try
    {
        CircuitTool__ oo = CircuitTool__(
            env, 
            pname, 
            path_input, 
            flag_random,
            flag_compute_output,
            flag_print_output,
            flag_tex,
            true
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