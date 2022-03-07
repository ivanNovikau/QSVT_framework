#ifndef QSP_H
#define QSP_H

#include "QSVT.h"

/**
 * @brief QSVT of the time evolution.
 * 
 */
class QSP__ : public QSVT__
{
    public:
        QSP__(const QuESTEnv& env, YCS fname_input, YCS path_inputs);
        ~QSP__();
        void read_angles();

    protected:
        void create_circuit_components();
        bool read_special_parameters(std::istringstream& iss, std::string& key_name);
        void save_basic_specific_data();
        void simulation();

    protected:
        /**
         * @brief Create the iterate W for the qubitization.
         */
        void create_W();

        /**
         * @brief Create the conjugate transposed iterator.
         */
        void create_iW();

    protected:
        YVQv t_grid_; // time grid
        uint32_t N_angles_;
        YVQv angles_phis_;
        YSCQ cu_; // controlled oracle circuit;
        YSCQ icu_; // controlled inversed oracle circuit;
        YSQ cw_;  // controlled W oracle;
        YSQ icw_; // inversed controlled W oracle;
};



#endif