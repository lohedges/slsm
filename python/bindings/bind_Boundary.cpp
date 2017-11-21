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

#include "Boundary.cpp"

using namespace slsm;

PYBIND11_MAKE_OPAQUE(std::vector<BoundaryPoint>)

void bind_Boundary(py::module &m)
{
    // STL containers.
    py::bind_vector<std::vector<BoundaryPoint>>(m, "VectorBoundaryPoint", py::module_local());

    // Class definition.
    py::class_<BoundaryPoint>(m, "BoundaryPoint", py::module_local(),
        "A data structure for boundary point data.")

        // Constructors.

        .def(py::init<>(), "Default constructor.")

        // Member data.

        .def_readonly("coord", &BoundaryPoint::coord,
            "Coordinate of the boundary point.")

        .def_readonly("normal", &BoundaryPoint::normal,
            "Inward pointing normal vector.")

        .def_readonly("length", &BoundaryPoint::length,
            "Integral length.")

        .def_readwrite("velocity", &BoundaryPoint::velocity,
            "Normal velocity.")

        .def_readwrite("negativeLimit", &BoundaryPoint::negativeLimit,
            "Movement limit in the negative direction.")

        .def_readwrite("positiveLimit", &BoundaryPoint::positiveLimit,
            "Movement limit in the positive direction.")

        .def_readwrite("isDomain", &BoundaryPoint::isDomain,
            "Whether the point lies within a grid spacing of the domain boundary.")

        .def_readwrite("isFixed", &BoundaryPoint::isFixed,
            "Whether the point is fixed.")

        .def_readonly("nSegments", &BoundaryPoint::nSegments,
            "The number of boundary segments that the point belongs to.")

        .def_readonly("segments", &BoundaryPoint::segments,
            "The indices of the two segments to which the point belongs.")

        .def_readonly("nNeighbours", &BoundaryPoint::nNeighbours,
            "The number of neighbouring boundary points.")

        .def_readonly("neighbours", &BoundaryPoint::neighbours,
            "The indices of the neighbouring points.")

        .def_readwrite("sensitivities", &BoundaryPoint::sensitivities,
            "The objective and constraint sensitivities.");

    // Class definition.
    py::class_<BoundarySegment>(m, "BoundarySegment", py::module_local(),
        "A data structure for boundary segment data.")

        // Constructors.

        .def(py::init<>(), "Default constructor.")

        // Member data.

        .def_readonly("start", &BoundarySegment::start,
            "Index of the start point.")

        .def_readonly("end", &BoundarySegment::end,
            "Index of the end point.")

        .def_readonly("element", &BoundarySegment::element,
            "The element cut by the boundary point.")

        .def_readonly("length", &BoundarySegment::length,
            "The length of the boundary segment.")

        .def_readonly("weight", &BoundarySegment::weight,
            "The weighting factor for the boundary segment.");

    // Class definition.
    py::class_<Boundary>(m, "Boundary", py::module_local(),
        "The discretised boundary of the level-set zero contour.")

        // Constructors.

        .def(py::init<>(), "Default constructor.")

        // Member functions.

        .def("discretise", &Boundary::discretise,
            "Use linear interpolation to compute the discretised boundary.",
            py::arg("levelSet"), py::arg("isTarget") = false)

        .def("computeNormalVectors", &Boundary::computeNormalVectors,
            "Compute the local normal vector at each boundary point.",
            py::arg("levelSet"))

        .def("computePerimeter", &Boundary::computePerimeter,
            "Compute the local perimeter for a boundary point.",
            py::arg("point"))

        // Member data.

        .def_readonly("points", &Boundary::points, "The vector of boundary points.")
        .def_readonly("segments", &Boundary::segments, "The vector of boundary segments.")
        .def_readonly("nPoints", &Boundary::nPoints, "The number of boundary points.")
        .def_readonly("nSegments", &Boundary::nSegments, "The number of boundary segments.")
        .def_readonly("length", &Boundary::length, "The total length of the boundary.");
}
