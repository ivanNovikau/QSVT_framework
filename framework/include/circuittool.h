#ifndef CIRCUITTOOL_H
#define CIRCUITTOOL_H

#include "BaseTool.h"


class CircuitTool__ : public BaseTool__
{
    public:
        /**
         * @param[in] pname  is a project name that defines names of input and output files.
         * @param[in] path_to_inputs is path to input files.
         * @param[in] flag_random create a random circuit;
         * @param[in] flag_compute_output if True, compute output states, otherwise compute only input states;
         * @param[in] flag_print_output if True, print output states on screen;
         * @param[in] flag_tex print or not .tex files;
         * @param[in] flag_hdf5  to create the hdf5 file;
         */
        CircuitTool__(
            const QuESTEnv& env, 
            YCS pname, 
            YCS path_to_inputs, 
            YCB flag_random,
            YCB flag_compute_output,
            YCB flag_print_output,
            YCB flag_tex,
            YCB flag_hdf5 = false
        );
        ~CircuitTool__();

        void launch();

    protected:
        void read_circuit_structure_from_file(YCS data);
        void create_random_circuit();

    protected:
        // to create a random circuit:
        bool flag_random_;
};

#endif