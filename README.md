# LibSLSM

<p>Copyright &copy; 2015-2017 <a href="http://lesterhedges.net">Lester Hedges</a>
<a href="http://www.gnu.org/licenses/gpl-3.0.html">
<img width="80" src="http://www.gnu.org/graphics/gplv3-127x51.png"></a></p>

[![Build Status](https://travis-ci.org/lohedges/slsm.svg?branch=master)](https://travis-ci.org/lohedges/slsm)

## About
LibSLSM is a small C++ library that forms an accompaniment to the paper
[Stochastic level-set method for shape optimisation](https://doi.org/10.1016/j.jcp.2017.07.010).
The library provides a simple set of building blocks for implementing
the method. Various example codes illustrate how to make make use of the
library and reproduce the data in the [paper](https://doi.org/10.1016/j.jcp.2017.07.010).

Unlike traditional level-set optimisation methods, which use a steepest
descent approach, SLSM updates the level set using a stochastic differential
equation. The addition of thermal noise enables the exploration of design
space during the optimisation process. In situations where the objective
function is non-convex, this allows for the possibility of escaping local
minima and converging to the true global optimum. See the
[paper](https://doi.org/10.1016/j.jcp.2017.07.010) for further details.

## Using the Code

LibSLSM is intended to be used on a [Unix-like](https://en.wikipedia.org/wiki/Unix-like)
operating system and has been extensively tested on macOS (formerly OS X) and
Linux. Windows users are advised to use the Ubuntu Bash shell on WIndows 10, or
[MinGW](http://www.mingw.org) to build LibSLSM. For modern, 64-bit Windows systems,
we recommend using [MYSYS2](http://msys2.github.io). This has been used to successfully
build LibSLSM on Windows 7 and 10. After installing MYSYS2 (following the instructions
on the [website](http://msys2.github.io)) you will need to install several additional
packages:

```bash
pacman -S cmake gcc git make python3
```

* Clone the LibSLSM source via:

```bash
git clone --recursive https://github.com/lohedges/slsm.git
```

(Alternatively, a zip file can be downloaded from
[here](https://github.com/lohedges/slsm/archive/master.zip).)

### Compiling and Installing:

* Change to the newly created `slsm` directory:

```bash
cd slsm
```

* Create a fresh build directory:

```bash
mkdir build
cd build
```
(This is recommended, although in source builds are allowed.)

* Configure, build, and install:

```bash
cmake .. && make -j4 install
```
(Here we use four cores for the build. If you have a lot of cores, try `make -j${nproc}`.)

On macOS you may wish to use LLVM's `libc++` rather than `libstdc++`. If so,
build as follows (this assumes a fresh build, or a clear cache):

```bash
CXXFLAGS=-stdlib=libc++ cmake .. && make -j4 install
```
By default, the library installs to `/usr/local`. Therefore, you may need admin
privileges for the final `make install` step above. An alternative is to change
the install location:

```bash
cmake -DCMAKE_INSTALL_PREFIX:PATH=MY_INSTALL_DIR .. && make -j4 install
```

(Note that there is no need to install the library in order to use it. You
can always build locally and link against the library using whatever path
is appropriate.)

### Linking with C/C++

To use LibSLSM with a C/C++ code first include the library header file
in the code:

```cpp
//example.cpp
#include <slsm/slsm.h>
```

Then to compile, we can use something like the following:
```bash
g++ -std=c++11 example.cpp -lslsm -lnlopt
```

This assumes that we have used the default install location `/usr/local/`. If
we specify an install location, we would use a command more like the following:

```bash
g++ -std=c++11 example.cpp -I/my/path/include -L/my/path/lib -lslsm -lnlopt
```

Note that the `-std=c++11` compiler flag is needed for `std::function` and `std::random`.

### External Dependencies
To aid portability, dependencies are handled via git
[submodlules](https://git-scm.com/book/en/v2/Git-Tools-Submodules). If you are using
LibSLSM for the first time, and didn't clone using the `--recursive` flag,
then you'll need to initialise the submodules:

```bash
git submodule update --init --recursive
```

Existing submodules can be updated using:

```bash
git submodule update --recursive
```

Brief details of the submodules are given below:

- The Optimise class makes use of [NLopt](http://ab-initio.mit.edu/wiki/index.php/NLopt).

- Python bindings are generated using [pybind11](https://github.com/pybind/pybind11).
See [Python](python/README.md) for details about using the Python module, pyslsm.

## Documentation

### Source code

Comprehensive documentation is provided via [Doxygen](www.doxygen.org). Make
sure you also have working [Graphviz](http://www.graphviz.org) (for
[DOT](http://www.graphviz.org/doc/info/lang.html)) and
[Tex Live](https://www.tug.org/texlive) (for [LaTeX](https://www.latex-project.org))
installations. To generate the documentation, run:

```bash
make doc
```

Following this, point your web browser to `doc/html/index.html`.

### Class Structure

For an overview of the class structure and details on how to instantiate
and use objects, see:
- [Class Overview](src/README.md)

### Tests

To learn how to compile and run unit tests, see:
- [Tests](tests/README.md)

### Examples
To get a feel for the how to write code using the library, see the
demonstration programs:
- [Demos](demos/README.md)

A set of utility scripts are provided for processing the output data:
- [Utils](utils/README.md)

### Python
Full Python bindings are generated using [pybind11](https://github.com/pybind/pybind11).
For more details, and to learn more about the Python extension module pyslsm, see:
- [Python](python/README.md)

## Acknowledgements
- Parts of this library were based on code written by
[Peter Dunning](http://www.abdn.ac.uk/engineering/people/profiles/peter.dunning).

- The Fast Marching Method implementation was adapted from
[Scikit-FMM](https://github.com/scikit-fmm/scikit-fmm).

## Disclaimer
Please be aware that this a working repository so the code should be used at
your own risk.

It would be great to hear from you if this code was of use in your research.
Email bugs, comments, and suggestions to lester.hedges+slsm@gmail.com.
