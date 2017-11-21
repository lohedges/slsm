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

namespace py = pybind11;

#include "Hole.cpp"

using namespace slsm;

void bind_Hole(py::module &m)
{
    // Class definition.
    py::class_<Coord>(m, "Coord", py::module_local(),
        "A data structure for two dimensional coordinates.")

        // Constructors.

        .def(py::init<>(), "Default constructor.")

        .def(py::init<double, double>(),
            "Constructor.", py::arg("x"), py::arg("y"))

        // Member data.

        .def_readwrite("x", &Coord::x, "The x position.")
        .def_readwrite("y", &Coord::y, "The y position.");

    // Class definition.
    py::class_<Hole>(m, "Hole", py::module_local(),
        "A data structure for circular holes.")

        // Constructors.

        .def(py::init<>(), "Default constructor.")

        .def(py::init<double, double, double>(),
            "Constructor", py::arg("x"), py::arg("y"), py::arg("r"))

        .def(py::init<const Coord&, double>(),
            " Constructor", py::arg("coord"), py::arg("r"))

        // Member data.

        .def_readwrite("coord", &Hole::coord, "The coordinates of the hole centre.")
        .def_readwrite("r", &Hole::r, "The radius of the hole.");

    // Class definition.
    py::class_<MutableFloat>(m, "MutableFloat", py::module_local(),
        "A mutable float to allow arguments to be passed by reference.")

		// Constructors.

        .def(py::init<>(), "Default constructor.")

        .def(py::init<double>(), "Constructor.", py::arg("value"))

        // Member data.

        .def_readwrite("value", &MutableFloat::value, "The value of the float.");
}
