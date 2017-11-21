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
#include <pybind11/stl_bind.h>

namespace py = pybind11;

#include "FastMarchingMethod.cpp"
#include "Heap.cpp"

using namespace slsm;

void bind_FastMarchingMethod(py::module &m)
{
    // Class definition.
    py::class_<FastMarchingMethod>(m, "FastMarchingMethod", py::module_local(),
        "Find approximate solitions to boundary value problems of the Eikonal equation.")

        // Constructors.

        .def(py::init<const Mesh&, bool>(),
            "Constructor.", py::arg("mesh"), py::arg("isTest") = false)

        // Member functions.

        .def("march", (void (FastMarchingMethod::*)(std::vector<double>&)) &FastMarchingMethod::march,
            "Reinitialise a signed distance function.",
            py::arg("signedDistance"))

        .def("march", (void (FastMarchingMethod::*)(std::vector<double>&,
            std::vector<double>&)) &FastMarchingMethod::march,
            "Extend boundary point velocities to nodes within the narrow band region.",
            py::arg("signedDistance"), py::arg("velocity"));
}
