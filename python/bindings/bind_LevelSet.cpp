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

#include "LevelSet.cpp"

using namespace slsm;

PYBIND11_MAKE_OPAQUE(std::vector<Coord>)
PYBIND11_MAKE_OPAQUE(std::vector<Hole>)

void bind_LevelSet(py::module &m)
{
    // STL containers.
    py::bind_vector<std::vector<Coord>>(m, "VectorCoord", py::module_local());
    py::bind_vector<std::vector<Hole>>(m, "VectorHole", py::module_local());

    // Class definition.
    py::class_<LevelSet>(m, "LevelSet", py::module_local(),
        "A two-dimensional level-set domain.")

        // Constructors.

        .def(py::init<unsigned int, unsigned int, double, unsigned int, bool>(),
            "Constructor.", py::arg("width"), py::arg("height"),
            py::arg("moveLimit") = 0.5, py::arg("bandWidth") = 6,
            py::arg("isFixedDomain") = false)

        .def(py::init<unsigned int, unsigned int, const std::vector<Hole>&, double,
            unsigned int, bool>(), "Constructor.", py::arg("width"),
            py::arg("height"), py::arg("holes"), py::arg("moveLimit") = 0.5,
            py::arg("bandWidth") = 6, py::arg("isFixedDomain") = false)

        .def(py::init<unsigned int, unsigned int, const std::vector<Coord>&, double,
            unsigned int, bool>(), "Constructor.", py::arg("width"),
            py::arg("height"), py::arg("points"), py::arg("moveLimit") = 0.5,
            py::arg("bandWidth") = 6, py::arg("isFixedDomain") = false)

        .def(py::init<unsigned int, unsigned int, const std::vector<Hole>&,
            const std::vector<Hole>&, double, unsigned int, bool>(),
            "Constructor.", py::arg("width"), py::arg("height"),
            py::arg("initialHoles"), py::arg("targetHoles"), py::arg("moveLimit") = 0.5,
            py::arg("bandWidth") = 6, py::arg("isFixedDomain") = false)

        .def(py::init<unsigned int, unsigned int, const std::vector<Hole>&,
            const std::vector<Coord>&, double, unsigned int, bool>(),
            "Constructor.", py::arg("width"), py::arg("height"),
            py::arg("initialHoles"), py::arg("targetPoints"), py::arg("moveLimit") = 0.5,
            py::arg("bandWidth") = 6, py::arg("isFixedDomain") = false)

        .def(py::init<unsigned int, unsigned int, const std::vector<Coord>&,
            const std::vector<Coord>&, double, unsigned int, bool>(),
            "Constructor.", py::arg("width"), py::arg("height"),
            py::arg("initialPoints"), py::arg("targetPoints"), py::arg("moveLimit") = 0.5,
            py::arg("bandWidth") = 6, py::arg("isFixedDomain") = false)

        // Member functions.

        .def("update", &LevelSet::update, "Update the level-set function."
            " The return value indicates whether the signed distance was reinitialised.",
            py::arg("timeStep"))

        .def("mask", (void (LevelSet::*)(const std::vector<Hole>&)) &LevelSet::mask,
            "Mask off a region of the domain.",
            py::arg("holes"))

        .def("mask", (void (LevelSet::*)(const std::vector<Coord>&)) &LevelSet::mask,
            "Mask off a region of the domain.",
            py::arg("points"))

        .def("reinitialise", &LevelSet::reinitialise,
            "Reinitialise the level set to a signed distance function.")

        .def("computeVelocities", (void (LevelSet::*)(const std::vector<BoundaryPoint>&))
            &LevelSet::computeVelocities,
            "Extend boundary point velocities to the level-set nodes.",
            py::arg("boundaryPoints"))

        .def("computeVelocities", (double (LevelSet::*)(std::vector<BoundaryPoint>&,
            MutableFloat&, const double, MersenneTwister&)) &LevelSet::computeVelocities,
            "Extend boundary point velocities to the level-set nodes."
            " Returns the time step scaling factor.",
            py::arg("boundaryPoints"), py::arg("timeStep"), py::arg("temperature"),
            py::arg("rng"))

        .def("computeGradients", &LevelSet::computeGradients,
            "Compute the modulus of the gradient of the signed distance function.")

        .def("computeAreaFractions", &LevelSet::computeAreaFractions,
            "Compute the material area fraction enclosed by the discretised boundary.")

        // Member variables.

        .def_readwrite("signedDistance", &LevelSet::signedDistance,
            "The nodal signed distance function.")

        .def_readwrite("velocity", &LevelSet::velocity,
            "The nodal normal velocity.")

        .def_readwrite("gradient", &LevelSet::gradient,
            "The nodal gradient of the signed distance function.")

        .def_readwrite("target", &LevelSet::target,
            "The target signed distance function.")

        .def_readonly("area", &LevelSet::area,
            "The total material area fraction.")

        .def_readonly("mesh", &LevelSet::mesh,
            "The fixed-grid mesh.")

        .def_readonly("moveLimit", &LevelSet::moveLimit,
            "The boundary movement limit (CFL condition).");
}
