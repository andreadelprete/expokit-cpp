
# expokit-cpp

The goal of this library is to provide a fast and comprehensive solution for computing the matrix exponential in C++. The original library [expokit](http://fortranwiki.org/fortran/show/Expokit) was written in Fortran.
Additionally is provided an integration utility using this type of algebraic calculation.

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


