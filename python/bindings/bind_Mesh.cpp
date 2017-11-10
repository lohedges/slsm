#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;

#include "Mesh.cpp"

using namespace slsm;

PYBIND11_MAKE_OPAQUE(std::vector<Element>)
PYBIND11_MAKE_OPAQUE(std::vector<Node>)

void bind_Mesh(py::module &m)
{
    // STL containers.
    py::bind_vector<std::vector<Element>>(m, "VectorElement", py::module_local());
    py::bind_vector<std::vector<Node>>(m, "VectorNode", py::module_local());

    // Class definition.
    py::class_<Element>(m, "Element", py::module_local(),
        "Data for an element in the two-dimensional fixed-grid mesh.")

        // Constructors

        .def(py::init<>(), "Constructor.")

        // Member data.

        .def_readonly("coord", &Element::coord,
            "The coordinates of the element centre.")

        .def_readonly("area", &Element::area,
            "The material area fraction of the element.");

    // Class definition.
    py::class_<Node>(m, "Node", py::module_local(),
        "Data for a node in the two-dimensional fixed-grid mesh.")

        // Constructors

        .def(py::init<>(), "Constructor.")

        // Member data.

        .def_readonly("coord", &Node::coord,
            "The coordinates of the node.")

        .def_readonly("neighbours", &Node::neighbours,
            "The indices of the neighbouring nodes.");

    // Class definition.
    py::class_<Mesh>(m, "Mesh", py::module_local(),
        "A two-dimensional fixed-grid mesh for the level-set domain.")

        // Constructors.

        .def(py::init<unsigned int, unsigned int>(),
            "Constructor.", py::arg("width"), py::arg("height"))

        // Member functions.

        .def("getClosestNode", (unsigned int (Mesh::*)(const Coord&) const) &Mesh::getClosestNode,
            "For a given coordinate, find the index of the closest node.",
            py::arg("coord"))

        .def("getClosestNode", (unsigned int (Mesh::*)(double, double) const) &Mesh::getClosestNode,
            "For a given coordinate, find the index of the closest node.",
            py::arg("x"), py::arg("y"))

        .def("getElement", (unsigned int (Mesh::*)(const Coord&) const) &Mesh::getElement,
            "For a given coordinate, find the element that contains that point.",
            py::arg("coord"))

        .def("getElement", (unsigned int (Mesh::*)(double, double) const) &Mesh::getElement,
            "For a given coordinate, find the element that contains that point.",
            py::arg("x"), py::arg("y"))

        // Member data.

        .def_readonly("nodes", &Mesh::nodes,
            "The nodes of the mesh.")

        .def_readonly("elements", &Mesh::elements,
            "The elements of the mesh.")

        .def_readonly("width", &Mesh::width,
            "The width of the mesh.")

        .def_readonly("height", &Mesh::height,
            "The height of the mesh.")

        .def_readonly("nElements", &Mesh::nElements,
            "The number of elements in the mesh.")

        .def_readonly("nNodes", &Mesh::nNodes,
            "The number of nodes in the mesh.");
}
