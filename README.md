# Diagonalisation of Quantum Observables (DoQO)

DoQO is a code for diagonalising obervables for quantum lattice models. It is
written in C/C++ and makes use of the SLEPc and PETSc libraries to run on
massively parallel distributed memory cluster. The code has been published along
with an accompanying paper in the journal of computer physics communications at
[http://www.sciencedirect.com/science/article/pii/S0010465511000099](http://www.sciencedirect.com/science/article/pii/S0010465511000099).

## To build DoQO. 

1. Ensure petsc and slepc are installed.

2. Create build folder. 
```
mkdir build
```

3. Run cmake to create makefiles. Specify boost location if not found automatically.
```
cd build
cmake ../
```
or 
```
cmake ../ -DBOOST_ROOT=<boost_path>
```
One can also specify an install prefix by setting the CMAKE_INSTALL_PREFIX variable with `-DCMAKE_INSTALL_PREFIX=<path>`

Not if one gets an error about "unset" from FindBoost, this is a cache issue, delete the build directory and start again.

4. Once cmake completes successfully build using.
```
make
```

5. Install binaries using
```
make install
```

