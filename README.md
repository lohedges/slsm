# LibLSM

<p>Copyright &copy; 2015-2016 <a href="http://lesterhedges.net">Lester Hedges</a>
<a href="http://www.gnu.org/licenses/gpl-3.0.html">
<img width="80" src="http://www.gnu.org/graphics/gplv3-127x51.png"></a></p>

## About
A simple C++ library to implement the Level Set Method (LSM) for peforming
structural topology optimisation. Based on code written by
[Peter Dunning](http://www.abdn.ac.uk/engineering/people/profiles/peter.dunning).

The aim is to provide a robust and exstensible, object-oriented topology
optimisation framework.

## Installation
A `Makefile` is included for building and installing LibLSM.

To compile LibLSM, then install the library, documentation, and demos:

```bash
$ make build
$ make install
```

By default, the library installs to `/usr/local`. Therefore, you may need admin
privileges for the final `make install` step above. An alternative is to change
the install location:

```bash
$ make PREFIX=MY_INSTALL_DIR install
```

Further details on using the Makefile can be found by running make without
a target, i.e.

```bash
$ make
```

Note that you don't have to install the library in order to use it. You can
always build locally then link against the library in the `lib` directory
inside of the repository.

## Compiling and linking
To use LibLSM with a C/C++ code first include the LibLSM header file somewhere
in the code.

```cpp
//example.cpp
#include <lsm/lsm.h>
```

Then to compile, we can use something like the following:

```bash
$ g++ -std=c++11 example.cpp -llsm -lnlopt
```

This assumes that we have used the default install location `/usr/local`. If
we specify an install location, we would use a command more like the following:

```bash
$ g++ -std=c++11 example.cpp -I/my/path/include -L/my/path/lib -llsm -lnlopt
```

Note that the `-std=c++11` compiler flag is needed for `std::function` and
`std::random`.

## Dependencies
Documentation can be generated using [Doxygen](http://www.stack.nl/~dimitri/doxygen)
by running `make doc`. Note that the `doc` target is a dependency of`install`,
so you will require a working doxygen installation. If you do not require
documentation then simply remove `doc` from the `install` dependency list in
the Makefile.

LibLSM uses the [Mersenne Twister](http://en.wikipedia.org/wiki/Mersenne_Twister)
psuedorandom number generator. A C++11 implementation using `std::random` is
included as a bundled header file, `MersenneTwister.h`. See the source code or
generate Doxygen documentation with `make doc` for details on how to use it.

The `Optimise` class makes use of [NLopt](http://ab-initio.mit.edu/wiki/index.php/NLopt).
Make sure that the library and header files are in your path. If not, do something like:

```bash
$ make OPTFLAGS="-I PATH_TO_NLOPT_HEADER -L PATH_TO_NLOPT_LIB" release
```

## Tests
A test suite is provided in the `tests` directory. To run all unit tests:

```bash
$ make test
```

Note that this will compile the tests (and library) using the default compilation
flags (`release`). To build the library and tests in a specific mode, run, e.g.

```bash
$ make devel test
```

## Demos
There are several example codes showing how to use the library:

* `demos/minimise_area.cpp`
* `demos/minimise_perimeter.cpp`
* `demos/minimise_area_constrained.cpp`
* `demos/minimise_perimeter_constrained.cpp`
* `demos/shape_match.cpp`
* `demos/dumbbell.cpp`
* `demos/umbrella_sample.cpp`

The Makefile will build standalone executables, which can be run (from the top level
directory) as follows, e.g.

```bash
$ ./demos/minimise_area

```

The demos serve as a test of our numerical implementation. The comments
in the source file provide further details of the calculations.

## Completed

### Mesh
Stores and initialises the two-dimensional finite element mesh. This
creates the required elements and nodes and stores information regarding their
connectivity. Support is provided for periodic and non-periodic meshes.

### LevelSet
Holds information relating to the level set function. Methods are
provided to initialise the signed distance function (either using a vector of
holes, or in a default "Swiss cheese" configuration) and to reinitialise it
using the fast marching method.

### FastMarchingMethod
An implementation of the fast marching method to find approximate solutions
to boundary value problems of the Eikonal equation. This is adapted from
[scikit-fmm](https://github.com/scikit-fmm/scikit-fmm), for which I've added
a few performance tweaks and bug fixes. The object provides functionality for
calculating signed distances and extension velocities.

### Boundary
An object for the discretised boundary. The `discretise` method solves for
a set of boundary points and segments given the current mesh and level set.
Following this, the `computeAreaFractions` method can determine the material
area fraction in each of the finite element cells.

### InputOutput
Provides functionality for reading and writing data structures. Currently
able to write: level set information, boundary points and segments, and
material area fractions. Data can be written as plain text (.txt) or in
ParaView readable VTK format. Methods should be able to read/write files
in the current directory, or from a user defined path.

### Optimise
Solve for the optimum boundary point velocity vector using the SLP method.
I have made numerous modifications to the method in order to remove unphysical
velocity capping and ensure that variable scaling is consistent and transparent.

### Sensitivity
A general purpose class for computing boundary point sensitivities by finite
differences. A callback function allows sensitivities to be calculated with
respect to an arbitrary function.

## Todo
* Work out why constraint fails when noisy sensitivity is used, e.g. perimeter
constraint. Minimising area at constant perimeter fails, whereas minimising
perimeter at constant area works.

## Limitations
* Limited to two-dimensional systems.
* Finite element mesh is assumed to be a fixed two-dimensional grid comprised
of square elements of unit side.
