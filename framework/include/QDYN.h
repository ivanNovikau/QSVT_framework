#ifndef QDYN_H
#define QDYN_H

#include "QSVT.h"

/**
 * @brief QSVT of the time evolution.
 * 
 */
class QDYN__ : public QSVT__
{
    public:
        QDYN__(const QuESTEnv& env, YCS fname_input, YCS path_inputs);
        void read_angles();

    protected:
        void create_circuit_components();
        bool read_special_parameters(std::istringstream& iss, std::string& key_name);
        void save_basic_specific_data();
        void simulation();

    protected:
        YVQv t_grid_; // time grid

};



#endif