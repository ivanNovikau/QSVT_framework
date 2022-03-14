#include "../include/BaseTool.h"
using namespace std;


BaseTool__::BaseTool__(
    const QuESTEnv& env, 
    YCS pname, 
    YCS path_to_inputs, 
    YCB flag_compute_output,
    YCB flag_circuit,
    YCB flag_tex,
    YCB flag_layers
) :
    env_(env),
    pname_(pname),
    path_inputs_(path_to_inputs),
    flag_compute_output_(flag_compute_output),
    flag_circuit_(flag_circuit),
    flag_tex_(flag_tex),
    flag_layers_(flag_layers)
{}


BaseTool__::~BaseTool__(){}


void BaseTool__::read_data()
{
    ifname_ = path_inputs_ + "/" + pname_ + format_file_;

    string current_path = filesystem::current_path();
    YMIX::print_log(env_, "Current path: " + current_path);
    YMIX::print_log(env_, "Start reading the input file: " + ifname_);

    string data;
    read_input_file(data);
    read_circuit_structure_from_file(data);
}


void BaseTool__::read_input_file(YS data)
{
    ifstream ff(ifname_);
    if(!ff.is_open()) throw "Error: Here there is not a file: " + ifname_;
    data = string((istreambuf_iterator<char>(ff)), istreambuf_iterator<char>());
    ff.close();

    // clean the buffer from empty lines and comments:
    istringstream istr(data);
    string data_clr = "";
    string line;
    while(getline(istr, line))
    {
        line = YMIX::remove_comment(line);
        line = YMIX::trim(line);
        if(line.find_first_not_of(' ') == string::npos)
            continue;
        data_clr += line + "\n";
    }
    std::transform(data_clr.begin(), data_clr.end(), data_clr.begin(), ::tolower);
    data = data_clr;
}