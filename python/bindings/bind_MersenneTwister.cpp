#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

namespace py = pybind11;

#include "MersenneTwister.h"

using namespace slsm;

void bind_MersenneTwister(py::module &m)
{
    // Class definition.
    py::class_<MersenneTwister>(m, "MersenneTwister", py::module_local(),
        "Mersenne-Twister pseudorandom number generator.")

        // Constructors.

        .def(py::init<>(), "Default constructor.")

        // Operator overloads.

		.def("__call__", &MersenneTwister::operator(),
            "Generate a uniform random number between 0 and 1.")

        // Member functions.

        .def("integer", &MersenneTwister::integer,
            "Generate and uniform random integer in the inclusive range.",
            py::arg("lower"), py::arg("upper"))

        .def("normal", (double (MersenneTwister::*)()) &MersenneTwister::normal,
            "Generate a normal distributed random number with zero mean and unit variance.")

        .def("normal", (double (MersenneTwister::*)(double, double)) &MersenneTwister::normal,
            "Generate a normal distributed random number with specified mean and variance).",
            py::arg("mean"), py::arg("variance"))

        .def("getSeed", &MersenneTwister::getSeed,
            "Get the value of the generator's seed.")

        .def("setSeed", &MersenneTwister::setSeed,
            "Set the value of the generator's seed.");
}
