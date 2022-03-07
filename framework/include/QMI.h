#ifndef QMI_H
#define QMI_H

#include "QSVT.h"

/**
 * @brief QSVT of the matrix inversion.
 * 
 */
class QMI__ : public QSVT__
{
    public:
        QMI__(const QuESTEnv& env, YCS fname_input, YCS path_inputs);
        void read_angles();

    protected:
        void create_circuit_components();
        bool read_special_parameters(std::istringstream& iss, std::string& key_name);
        void save_basic_specific_data();
        void simulation();
};



#endif