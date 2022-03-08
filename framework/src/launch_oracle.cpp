#include "../include/oracletool.h"
using namespace std;

/**
 * @brief To launch the oracletool.
 * @param argv oracletool, [project_name], [path_to_input_files], [flag_compute_output]
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

    // flag to compute output from an oracle
    bool flag_compute_output = true;
    if(argc > 3) flag_compute_output = stoi(string (argv[3]));
    if(flag_compute_output) YMIX::print_log(env, "Flag compute output oracle = true");
    else                    YMIX::print_log(env, "Flag compute output oracle = false");
    
    try
    {
        // create the tool for the analysis of an oracle
        OracleTool__ oo = OracleTool__(
            env, 
            pname, 
            path_input, 
            flag_compute_output
        );
        
        // launch the circuit with different input states:
        oo.launch();
    }
    catch(const string& e)
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