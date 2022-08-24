# QSVT_framework
QSVT and QSP Hamiltonian simulations and matrix inversion.
If one wants to check the code used in Ref. [1], go to the branch PRA22.

## Prepare the code
One needs `cmake` to compile the framework and `python` to analyse the simulations.
In particular, one needs the `qiskit` and `h5py` modules:

`pip install qiskit h5py`

Remark: the code does not support multi-GPU parallelization.

## Compilation
The framework uses the QuEST code (https://github.com/quest-kit/QuEST) as a submodule. Therefore, while cloning, use the flag `--recurse-submodules`.
Since the code also needs the HDF5 library, in the `./framework/external/QuEST/CMakeLists.txt` file, add the line `find_package(HDF5 REQUIRED COMPONENTS CXX HL)` and replace the line `target_link_libraries(${OUTPUT_EXE} QuEST m)` by the line `target_link_libraries(${OUTPUT_EXE} QuEST m HDF5::HDF5)`.

Correct the `-DGPU_COMPUTE_CAPABILITY` in the `tobuild` files according to your GPU device
(https://developer.nvidia.com/cuda-gpus#compute).

To compile the oracle-tool:
1. Go to the folder `./framework/build_oracle`.
2. `source tobuild`
3. `make`

To compile the framework:
1. Go to the folder `./framework/build_framework`.
2. `source tobuild`
3. `make`

## Run the oracle-tool
The oracle-tool, `oracletool`, reads the `.oracle` file to constrcut the circuit described there, calculates the output states of this circuit.
If necessary, `oracletool` builds the `.circuit` and `.tex` representation of the circuit.

To run the oracle-tool, use the following command:

`[path_to_QSVT_framework]/framework/build_oracle/oracletool [oracle-name] [work-directory] [optional-parameters]`

The `oracletool` searches for the file `[oracle-name].oracle` in the folder `[work-directory]`.

The optional parameters, `[optional-parameters]`, include : 

`-flag_compute_output`: if `-flag_compute_output 1`, then `oracletool` calculates output states for the specified input states from the `.oracle` file.
By default, `-flag_compute_output 1`.

`-flag_print_output`: to print output states on screen.
By default, `-flag_print_output = 1`.

`-flag_circuit`: to write the `.circuit` files. 
By default, `-flag_circuit 1`.

`-flag_tex`: to write the `.tex` files. 
By default, `-flag_tex 1`.

`-tex_CL`: length of each row in the circuit, when it is printed to the `.tex` file.
By default, `-tex_CL 10`.

`-flag_print_zero_anc`: print only states, where all ancillae are in the zero state
(assuming that all ancilla qubit are the high-priority qubits).

To understand the format of the `.oracle` file, see the wiki page:<br> 
https://github.com/ivanNovikau/QSVT_framework/wiki

## Run the circuit reader:
The circuit reader, `qc_circuit`,<br> 
--> either reads a `.circuit` file to create a quantum circuit and to compute its output states (`qc_circuit` always takes the zero input state);<br> 
--> or creates a random circuit, computes its output states and creates the `.circuit` output file of this random circuit.

To run the circuit reader, use the following command:

`[path_to_QSVT_framework]/framework/build_circuit/qc_circuit [file-name] [work-directory] [optional-parameters]`

`[optional-parameters]` include the same parameters as the `oracletool` (except the `-flag_circuit`).<br>
`[optional-parameters]` also include `-flag_random`: <br>
      if `-flag_random 0`, then `qc_circuit` reads the `[file-name].circuit` file;<br>
      if `-flag_random 1`, then `qc_circuit` reads the `[file-name].random` file. The `[file-name].random` file contains the number of qubits and the number of gates in the random circuit to create. Then, `qc_circuit` creates the random circuit and writes it down to the `[file-name].circuit` file.


## Run the framework
The framework, `qsvt`, takes the `.oracle` files of the block-encoded matrix and of the initialization circuit and produces the QSVT (QSP) circuit for the parameters described in the `.qsp` file and using the rotation angles from the `.angles` file(s).

To run the framework, use the following command

`[path_to_QSVT_framework]/framework/qsvt [oracle-name] [work-folder with the .qsp file] [case-to-sim]`

The framework searches for the following files in the directory `[work-folder with the .qsp file]`:

`[oracle-name].qsp`: contains information about the QSVT(QSP) parameters such as the time interval for the Hamiltonian simulation;

`[oracle-name].oracle`: contains the circuit, which encodes the Hamiltonian of the simulated system;

The QSVT (QSP) circuit is initialized with a quantum state, which should be read from the file
`[oracle-name].init_state` if `sel_init vector [matrix-norm]` is set in the `.qsp` file.
If `sel_init oracle [matrix-norm]`, then the circuit from the `[oracle-name]_init.oracle` file is used to initialize the QSVT circuit.

`.angles`: file(s), which contain(s) the rotation angles for the QSP (QSVT) approximation.

The `.oracle` file might also need one or several file(s) `.condR_profile`, which contain profiles of angles for the condition rotation gates.

The framework can simulate QSP Hamiltonian simulations (only for hermitian Hamiltonians), QSVT Hamiltonian simulation, QSVT matrix inversion:<br>
`[case-to-sim] = qsp`: QSP Hamiltonian simulation;<br>
`[case-to-sim] = qsvt-dyn`: QSVT Hamiltonian simulation;<br>
`[case-to-sim] = qsvt-mi`: QSVT matrix inversion.

For `[case-to-sim] = qsp`, the `.angles` file has the following name:

`angles_t[t*1000]_eps[-log10(eps)].angles`,

where `t` is the normalized time interval to simulate (`t = t_orig*[matrix-norm]`, where `[matrix-norm] = H/H_norm`: `|H_norm| <= 1`).

For `[case-to-sim] = qsvt-dyn`, one needs two `.angles` files:

`angles_even_param[t*1000]_eps[-log10(eps)].angles` (to approximate `cos`)<br>
`angles_odd_param[t*1000]_eps[-log10(eps)].angles`  (to approximate `i*sin`).

For `[case-to-sim] = qsvt-mi`, the `.angles` file has the following name:

`angles_odd_param[kappa*1000]_eps[-log10(eps)].angles`,

where `kappa` is the condition number of the matrix to invert.
The condition number is found as the ratio between the max. and min. singular values.

For `qsp`, the angles can be calculated by the following code: https://github.com/microsoft/Quantum-NC/tree/main/src/simulation/qsp

For `qsvt`, the angles can be calculated by the following code: https://github.com/qsppack/QSPPACK

The framework does not create the `.tex` description of the corresponding circuits and subcircuits, but only some `.circuit` files if `flag_circuit true` in the `.qsp` file. By default, `flag_circuit false`.
It is highly recommended to keep `flag_circuit false` since the size of the resulting `.circuit` files for a whole QSVT (QSP) circuit might be of dozens of GB, and the writing of the `.circuit` file can slow down the calculations.

To produce the `.tex` files for the QSVT circuit, one needs to create the `.circuit` files and then use the `qc_circuit` program to create the `.tex` files.

Remark: the framework does not produce the `.circuit` files of the block-encoded matrix (described in the `.oracle` files) even if `flag_circuit true`.
To create the `.circuit` and `.tex` files of the `.oracle` circuits, one should use the `oracletool`.


# Simulations

## Cold linear plasma XW waves in 1-D space
The corresponding simulations using the QSP with a hermitian Hamiltonian are given in the following directories:

`./simulations/XW_waves/small-case`
`./simulations/XW_waves/ref-case`

The `small-case` is the simulation with `Nx = 64` spatial points (`nj 6` in `small-case\*.oracle` files).
This simulation can be used to test the QSP computation.

The `ref-case`, which is presented in Ref. 1, has `Nx = 1024` spatial points (`nj 10` in `ref-case\*.oracle` files).




# References
1. QSP for simulating cold plasma waves: https://arxiv.org/abs/2112.06086<br>
2. QuEST quantum simulator: https://www.nature.com/articles/s41598-019-47174-9<br>
3. Hamiltonian Simulation by Qubitization: https://quantum-journal.org/papers/q-2019-07-12-163/<br>
4. Quantum Singular Value Transformation (QSVT): https://dl.acm.org/doi/10.1145/3313276.3316366<br>
5. Calculation of the QSP angles: https://quantum-journal.org/papers/q-2019-10-07-190/<br>
6. Calculation of the QSVT angles: https://journals.aps.org/pra/abstract/10.1103/PhysRevA.103.042419<br>


























