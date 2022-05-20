
# expokit-cpp

This library provides fast algorithms for computing the exponential of dense matrices in C++.

Class [MatrixExponential](https://github.com/andreadelprete/expokit-cpp/blob/master/include/MatrixExponential.hpp):
* an optimized version of the scaling-and-squaring algorithm (Higham 2005 "The scaling and squaring method for the matrix exponential revisited") for computing the matrix exponential based on the code found in [Eigen3-unsupported](https://eigen.tuxfamily.org/dox/unsupported/group__MatrixFunctions__Module.html#matrixbase_exp);
* options for the user to customize the amount of computation performed by the algorithm to achieve the desired performance-accuracy trade-off (e.g. specifying the maximum number of matrix-matrix multiplications);
* an algorithm for computing the product between the matrix exponential and a given vector/matrix, which is more efficient than computing the matrix exponential and multiplying it times the vector/matrix (assuming the matrix is skinny);
* a balancing algorithm to speed-up the matrix exponential computation;

Class [LDSUtility](https://github.com/andreadelprete/expokit-cpp/blob/master/include/LDSUtility.hpp):
* a utility class for computing the evolution of the state of a linear dynamical system, and its first two integrals, all of which are based on the computation of matrix exponentials.

All the C files:
* a C version of the library [expokit](http://fortranwiki.org/fortran/show/Expokit) (Sidje, 1998 "Expokit: A Software Package for Computing Matrix Exponentials", ACM Transactions on Mathematical Software), which is a Fortran library that we translated to C using `f2c`.

In [this file](https://github.com/andreadelprete/expokit-cpp/blob/master/tests/test_small.cpp) you can find a test comparing all the implemented methods for computing the matrix exponential. We have found that for the matrix sizes we were interested in, our optimized version of the scaling-and-squaring algorithm based on [Eigen3-unsupported](https://eigen.tuxfamily.org/dox/unsupported/group__MatrixFunctions__Module.html#matrixbase_exp) was the fastest. For example, using a matrix size of 24, these were the computation times we obtained:
```
Our MatrixExponential class: 0.068 ms
Eigen (from eigen3-unsupported): 0.079 ms
DGPADM (from expokit): 0.160 ms
DGCHBV (from expokit): 0.094 ms

Our MatrixExponential class - computeExpTimesVector: 0.044 ms
```

## Dependencies
1. `f2c` - used in porting routines from Fortran
2. `eigen3` - Linear Algebra library
3. `lapack` - Linear Algebra library

## Build

First of all install `f2c` :

`sudo apt install f2c`

Then install `eigen3`, which can be installed from different sources. We propose here to install it from [robotpkg](http://robotpkg.openrobots.org/debian.html) using:

`sudo apt install robotpkg-eigen3`

To compile expokit you need to provide to `pkg-conifg` the path both to `eigen3.pc` and `f2c.pc`. The location of the former depends on how `eigen3` was installed. We are assuming here it was installed from `robotpkg`. 

Regarding `f2c`, the current package provided by `apt` does not include the `f2c.pc` file, but we provide our own that can be found in this repo at `pkg-config/f2c.pc`. 

The `lapack` dependency should be resolved by the build system. Summing up your `PKG_CONFIG_PATH` should look like this:
```
export PKG_CONFIG_PATH=/opt/openrobots/lib/pkgconfig:$PKG_CONFIG_PATH
export PKG_CONFIG_PATH=/opt/openrobots/share/pkgconfig:$PKG_CONFIG_PATH
export PKG_CONFIG_PATH=<path to f2c.pc>:$PKG_CONFIG_PATH
```
Substitute `<path to f2c.pc>` with your path, it can be located anywhere.

## Citing Expokit-cpp

To cite Expokit-cpp in your academic research, please use the following bibtex line:

```
@article{Hammoud2022,
author = {Hammoud, Bilal and Olivieri, Luca and Righetti, Ludovic and Carpentier, Justin and {Del Prete}, Andrea},
doi = {10.1007/s11044-022-09818-z},
issn = {1573-272X},
journal = {Multibody System Dynamics},
number = {4},
pages = {443--460},
title = {{Exponential integration for efficient and accurate multibody simulation with stiff viscoelastic contacts}},
url = {https://doi.org/10.1007/s11044-022-09818-z},
volume = {54},
year = {2022}
}
```


