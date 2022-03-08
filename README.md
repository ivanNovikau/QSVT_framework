# QSVT_framework
QSVT and QSP Hamiltonian simulations and matrix inversion

## Prepare the code
The code uses QuEST quantum emulator (https://github.com/quest-kit/QuEST).
One needs to clone the QuEST repository.
Then, create the global bash variable $QuESTHome to store the absolute path to the QuEST/ folder.

One also needs `cmake` to compile the framework and `python` to analyse the simulations.

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

## Run the framework
To run the oracle-tool, use the following command:

`[path_to_QSVT_framework]/framework/build_oracle/oracletool [name of oracle] [path to the input file] [flag_output]`

The `oracletool` searches for the file `[name of oracle].oracle` in the folder `[path to the input file]`.

The flag `[flag_output] = 1` - then calculate output states for the specified input states from the `.oracle` file;<br> 
`[flag_output] = 0` - then do not compute the output states (just produce `.circuit` and `.tex` files).


