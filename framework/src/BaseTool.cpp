#include "../include/BaseTool.h"
using namespace std;


BaseTool__::BaseTool__(
    const QuESTEnv& env, 
    YCS pname, 
    YCS path_to_inputs, 
    YCB flag_compute_output,
    YCB flag_print_output,
    YCB flag_circuit,
    YCB flag_tex,
    YCB flag_layers,
    YCB flag_hdf5
) :
    env_(env),
    pname_(pname),
    path_inputs_(path_to_inputs),
    flag_compute_output_(flag_compute_output),
    flag_print_output_(flag_print_output),
    flag_circuit_(flag_circuit),
    flag_tex_(flag_tex),
    flag_layers_(flag_layers),
    flag_hdf5_(flag_hdf5)
{
    string current_path = filesystem::current_path();
    YMIX::print_log("Current path: " + current_path);

    if(flag_hdf5_)
    {
        string filename_out = path_inputs_ + "/" + pname_ + ENDING_FORMAT_OUTPUT;

        YMIX::print_log("Creating the output .hdf5 file:\n\t "s + filename_out);
        hfo_.create(filename_out);
        hfo_.add_group("basic"); 
        hfo_.add_group("states"); 
        hfo_.add_group("constants"); 

        // date of simulation:
        string str_date_time;
        YMIX::get_current_date_time(str_date_time);
        hfo_.add_scalar(str_date_time, "date-of-simulation", "basic");

        // save the project name and the path to input files:
        hfo_.add_scalar(pname_, "project-name", "basic");
        hfo_.add_scalar(path_inputs_, "path-inputs", "basic");
        hfo_.add_scalar(filesystem::current_path(), "launch-path", "basic");

        // close the file:
        hfo_.close();
    }
}


BaseTool__::~BaseTool__(){}


void BaseTool__::read_data()
{
    ifname_ = path_inputs_ + "/" + pname_ + format_file_;
    YMIX::print_log("\n----------------------------------------------");
    YMIX::print_log("Start reading the input file: " + ifname_);

    try
    {
        string data;
        read_input_file(data);
        read_circuit_structure_from_file(data);
    }
    catch(YCS e)
    {
        throw "Error while reading the file["s + ifname_ + "]:\n"s + e;
    }
}


void BaseTool__::read_input_file(YS data, YCS file_name)
{   
    string file_name_res = file_name;
    if(file_name_res.empty()) file_name_res = ifname_;

    ifstream ff(file_name_res);
    if(!ff.is_open()) throw "Error: Here there is not a file: " + file_name_res;
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
    // std::transform(data_clr.begin(), data_clr.end(), data_clr.begin(), ::tolower);
    data = data_clr;
}