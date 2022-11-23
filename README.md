# QSVT_framework
QSVT Hamiltonian simulations and matrix inversion.
To check the code used in Ref. [1](https://github.com/ivanNovikau/QSVT_framework/wiki/References), go to the branch PRA22.

## Prepare the code
One needs `cmake` to compile the framework and `python` to analyse the simulations.
In particular, one needs the `qiskit` and `h5py` modules:

`pip install qiskit h5py`

Remark: the code does not support multi-GPU parallelization.

## Compilation
The framework uses the [QuEST code](https://github.com/quest-kit/QuEST), Ref. [2](https://github.com/ivanNovikau/QSVT_framework/wiki/References) as a submodule. Therefore, while cloning, use the flag `--recurse-submodules`.
Since the code also needs the HDF5 library, in the `./framework/external/QuEST/CMakeLists.txt` file, add the line `find_package(HDF5 REQUIRED COMPONENTS CXX HL)` and replace the line `target_link_libraries(${OUTPUT_EXE} QuEST m)` by the line `target_link_libraries(${OUTPUT_EXE} QuEST m HDF5::HDF5)`.

Correct the `-DGPU_COMPUTE_CAPABILITY` in the `tobuild` files according to your GPU device
(https://developer.nvidia.com/cuda-gpus#compute).

To compile the oracle-tool:
1. Go to the folder `./framework/build_oracle`.
2. `source tobuild`
3. `make`

## Run the oracle-tool
The oracle-tool, `oracletool`, reads the `.oracle` file to constrcut the circuit described there, calculates the output states of this circuit.
If necessary, `oracletool` builds the `.circuit` and `.tex` representation of the circuit.

To run the oracle-tool, use the following command:

`[path_to_QSVT_framework]/framework/build_oracle/oracletool [oracle-name] [work-directory]`

The `oracletool` searches for the file `[oracle-name].oracle` in the folder `[work-directory]`.

To understand the format of the `.oracle` file, see the [wiki](https://github.com/ivanNovikau/QSVT_framework/wiki) page. 


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




























