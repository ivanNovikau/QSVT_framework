#include "../include/oracletool.h"
using namespace std;

/**
 * @brief To launch the oracletool.
 * @param argv oracletool [project_name] [path_to_input_files]
 */
int main(int argc, char *argv[])
{
    uint32_t id_arg;
    QuESTEnv env = createQuESTEnv();
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

    // --- project name ---
    id_arg = 1;
    string pname(argv[id_arg]);
    
    // --- path to input files ---
    id_arg += 1;
    string path_input(argv[id_arg]);

    // file of the log file:
    YMIX::LogFile::name_global_ = path_input + "/" + pname + ".clog";
    YMIX::LogFile cf(true);

    cout << "\n\n";
    cout << "Project name: "        << pname << endl;
    cout << "Path to input files: " << path_input << "/" << endl;

    try
    {
        OracleTool__ oo = OracleTool__(
            env, 
            pname, 
            path_input
        );
        oo.launch();
    }
    catch(YCS e)
    {
        std::cerr << "\n" << e << endl;
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