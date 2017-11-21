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

PYBIND11_MAKE_OPAQUE(std::vector<int>)
PYBIND11_MAKE_OPAQUE(std::vector<unsigned int>)
PYBIND11_MAKE_OPAQUE(std::vector<double>)
PYBIND11_MAKE_OPAQUE(std::vector<bool>)

void bind_Boundary(py::module &);
void bind_FastMarchingMethod(py::module &);
void bind_Hole(py::module &);
void bind_InputOutput(py::module &);
void bind_LevelSet(py::module &);
void bind_MersenneTwister(py::module &);
void bind_Mesh(py::module &);
void bind_Optimise(py::module &);
void bind_Sensitivity(py::module &);

PYBIND11_MODULE(pyslsm, m)
{
    // STL containers.
    py::bind_vector<std::vector<int>>(m, "VectorInt", py::module_local());
    py::bind_vector<std::vector<unsigned int>>(m, "VectorUnsignedInt", py::module_local());
    py::bind_vector<std::vector<double>>(m, "VectorDouble", py::module_local());
    py::bind_vector<std::vector<bool>>(m, "VectorBool", py::module_local());

    // Class bindings.
    bind_Boundary(m);
    bind_FastMarchingMethod(m);
    bind_Hole(m);
    bind_InputOutput(m);
    bind_LevelSet(m);
    bind_MersenneTwister(m);
    bind_Mesh(m);
    bind_Optimise(m);
    bind_Sensitivity(m);
}
