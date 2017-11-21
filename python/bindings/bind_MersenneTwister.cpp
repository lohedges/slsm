/*
  Copyright (c) 2015-2017 Lester Hedges <lester.hedges+slsm@gmail.com>

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

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
