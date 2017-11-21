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

#include "Optimise.cpp"

using namespace slsm;

void bind_Optimise(py::module &m)
{
    // Class definition.
    py::class_<Optimise>(m, "Optimise", py::module_local(),
        "Solve for the optimum velocity vector.")

        // Constructors.

        .def(py::init<std::vector<BoundaryPoint>&, std::vector<double>&, std::vector<double>&,
            MutableFloat&, double, bool, const std::vector<bool>&>(), "Constructor.",
            py::arg("boundaryPoints"), py::arg("constraindDistances"), py::arg("lambdas"),
            py::arg("timeStep"), py::arg("maxDisplacement") = 0.5, py::arg("isMax") = false,
            py::arg("isEquality") = std::vector<bool>())

        // Member functions.

        .def("solve", &Optimise::solve,
            "Execute the NLopt solver to find the optimium velocity vector."
            " Returns the optimum change in the objective function.")

        .def("queryReturnCode", &Optimise::queryReturnCode,
            "Query the NLopt return code.");
}
