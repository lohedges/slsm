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
