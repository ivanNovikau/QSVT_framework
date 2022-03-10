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
To run the oracle-tool, use the following command:

`[path_to_QSVT_framework]/framework/build_oracle/oracletool [name of oracle] [path to the input file] [flag_output]`

The `oracletool` searches for the file `[name of oracle].oracle` in the folder `[path to the input file]`.

`[flag_output] = 1`: calculate output states for the specified input states from the `.oracle` file;<br> 
`[flag_output] = 0`: do not compute the output states (just produce `.circuit` and `.tex` files).

To understand the format of the `.oracle` file, see the wiki page:<br> 
https://github.com/ivanNovikau/QSVT_framework/wiki

## Run the framework
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

## References
1. QSP for simulating cold plasma waves: https://arxiv.org/abs/2112.06086<br>
2. QuEST quantum simulator: https://www.nature.com/articles/s41598-019-47174-9<br>
3. Hamiltonian Simulation by Qubitization: https://quantum-journal.org/papers/q-2019-07-12-163/<br>
4. Quantum Singular Value Transformation (QSVT): https://dl.acm.org/doi/10.1145/3313276.3316366<br>
5. Calculation of the QSP angles: https://quantum-journal.org/papers/q-2019-10-07-190/<br>
6. Calculation of the QSVT angles: https://journals.aps.org/pra/abstract/10.1103/PhysRevA.103.042419<br>























