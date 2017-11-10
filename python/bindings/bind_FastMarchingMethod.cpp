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
