# QSVT_framework
QSVT and QSP Hamiltonian simulations and matrix inversion

## Prepare the code
The code uses QuEST quantum emulator (https://github.com/quest-kit/QuEST).
One needs to clone the QuEST repository.
Then, create the global bash variable $QuESTHome to store the absolute path to the QuEST/ folder.

One also needs `cmake` to compile the framework and `python` to analyse the simulations.

## Compilation
To compile the oracle-tool:
1. Go to the folder `./framework/build_oracle`.
2. `source tobuild`
3. `make`

To compile the framework:
1. Go to the folder `./framework/build_framework`.
2. `source tobuild`
3. `make`
