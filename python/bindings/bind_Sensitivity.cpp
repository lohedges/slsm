#include <pybind11/pybind11.h>
#include <pybind11/functional.h>

namespace py = pybind11;

#include "Sensitivity.cpp"

using namespace slsm;

void bind_Sensitivity(py::module &m)
{
    // Class definition.
    py::class_<Callback>(m, "Callback", py::module_local(),
        "Calculate the value of the function for a small displacement of the boundary point.")

        // Constructors.

        .def(py::init<>(), "Default constructor.")

        // Member data.

        .def_readwrite("callback", &Callback::callback, "The callback function.");

    // Class definition.
    py::class_<Sensitivity>(m, "Sensitivity", py::module_local(),
        "Calculate finite-difference boundary point sensitivities.")

        // Constructors.

        .def(py::init<double>(), "Constructor.", py::arg("delta") = 1e-4)

        // Member functions.

        .def("computeSensitivity", &Sensitivity::computeSensitivity,
            "Compute the finite-difference sensitivity for an arbitrary function.",
            py::arg("point"), py::arg("callback"))

        .def("itoCorrection", (void (Sensitivity::*)
            (Boundary&, const LevelSet&, double) const) &Sensitivity::itoCorrection,
            "Apply deterministic Ito correction to the objective sensitivity.",
            py::arg("boundary"), py::arg("levelSet"), py::arg("temperature"))

        .def("itoCorrection", (void (Sensitivity::*)
            (Boundary&, double) const) &Sensitivity::itoCorrection,
            "Apply deterministic Ito correction to the objective sensitivity.",
            py::arg("boundary"), py::arg("temperature"));
}
