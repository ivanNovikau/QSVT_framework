# QSVT_framework
QSVT and QSP Hamiltonian simulations and matrix inversion

## Prepare the code
The code uses QuEST quantum emulator (https://github.com/quest-kit/QuEST).
One needs to clone the QuEST repository.
Then, create the global bash variable $QuESTHome to store the absolute path to the QuEST/ folder.

One also needs `cmake` to compile the framework and `python` to analyse the simulations.
In particular, one needs the `qiskit` and `h5py` modules:

`pip install qiskit h5py`

Correct the `-DGPU_COMPUTE_CAPABILITY` in the `tobuild` files according to your GPU device
(https://developer.nvidia.com/cuda-gpus#compute).

Remark: the code does not support multi-GPU parallelization.

## Compilation
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

`[path_to_QSVT_framework]/framework/build_oracle/oracletool [name of oracle] [path to the input file] [flag_output] [flag-circuit] [flag-tex] [tex-circuit-length]`

The `oracletool` searches for the file `[name of oracle].oracle` in the folder `[path to the input file]`.

`[flag_output]`, optional: if `1`, then calculate output states for the specified input states from the `.oracle` file;<br> 
if `0`, then do not compute the output states (just produce `.circuit` and `.tex` files).<br>
By default, `[flag_output] = 1`;

`[flag-circuit]`, optional: do print or not the `.circuit` files. 
By default, `[flag-circuit] = 1`.

`[flag-tex]`, optional: do print or not the `.tex` files. 
By default, `[flag-tex] = 1`.

`[tex-circuit-length]`, optional: length of each row in the circuit, when it is printed to the `.tex` file (by default, ).<br>
By default, `[tex-circuit-length] = 10`.

To understand the format of the `.oracle` file, see the wiki page:<br> 
https://github.com/ivanNovikau/QSVT_framework/wiki

## Run the circuit constructor:
The circuit constructor, `qc_circuit`, can read one or several `.circuit` files to combine them into a single quantum circuit and calculate the corresponding output state (`qc_circuit` always takes the zero input state).
The `qc_circuit` connects the circuits by comparing the register names.
The program can produce the `.circuit` and `.tex` representations of the resulting circuit.

## Run the framework
The framework, `qsvt`, takes the `.oracle` files of the block-encoded matrix and of the initialization circuit and produces the QSVT (QSP) circuit for the parameters described in the `.qsp` file and using the rotation angles from the `.angles` file(s).

To run the framework, use the following command

`[path_to_QSVT_framework]/framework/qsvt [name of oracle] [work-folder with the .qsp file] [case-to-sim]`

The framework searches for the following files in the directory `[work-folder with the .qsp file]`:

`[name of oracle].qsp`: contains information about the QSVT(QSP) parameters such as the time interval for the Hamiltonian simulation;

`[name of oracle].oracle`: contains the circuit, which encodes the Hamiltonian of the simulated system;

`[name of oracle].init_state`: contains the initial state vector;

`.angles`: file(s), which contain(s) the rotation angles for the QSP (QSVT) approximation.

The `[name of oracle].oracle` file might also need one or several file(s) `.condR_profile`, which contains a profile of angles for the condition rotation gates.

The framework can simulate QSP Hamiltonian simulations (only for hermitian Hamiltonians), QSVT Hamiltonian simulation, QSVT matrix inversion:<br>
`[case-to-sim] = qsp`: QSP Hamiltonian simulation;<br>
`[case-to-sim] = qsvt-dyn`: QSVT Hamiltonian simulation;<br>
`[case-to-sim] = qsvt-mi`: QSVT matrix inversion.

For `[case-to-sim] = qsp`, the `.angles` file has the following name:

`angles_t[t*1000]_eps[-log10(eps)].angles`,

where `t` is the normalized time interval to simulate (`t = t_orig*norm`, where `norm = H/H_norm`: `|H_norm| <= 1`).

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

## References
1. QSP for simulating cold plasma waves: https://arxiv.org/abs/2112.06086<br>
2. QuEST quantum simulator: https://www.nature.com/articles/s41598-019-47174-9<br>
3. Hamiltonian Simulation by Qubitization: https://quantum-journal.org/papers/q-2019-07-12-163/<br>
4. Quantum Singular Value Transformation (QSVT): https://dl.acm.org/doi/10.1145/3313276.3316366<br>
5. Calculation of the QSP angles: https://quantum-journal.org/papers/q-2019-10-07-190/<br>
6. Calculation of the QSVT angles: https://journals.aps.org/pra/abstract/10.1103/PhysRevA.103.042419<br>























