#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "InputOutput.cpp"

using namespace slsm;

void bind_InputOutput(py::module &m)
{
    // Class definition.
    py::class_<InputOutput>(m, "InputOutput", py::module_local(),
        "Functionality for reading and writing level set data.")

        // Constructors.

        .def(py::init<>(), "Constructor.")

        // Member functions.

        .def("saveLevelSetVTK", (void (InputOutput::*)(const unsigned int&,
            const LevelSet&, bool, bool, const std::string&) const) &InputOutput::saveLevelSetVTK,
            "Write the level set to a ParaView VTK file.",
            py::arg("datapoint"), py::arg("levelSet"), py::arg("isVelocity") = false,
            py::arg("isGradient") = false, py::arg("outputDirectory") = "")

        .def("saveLevelSetVTK", (void (InputOutput::*)(const std::string&,
            const LevelSet&, bool, bool) const) &InputOutput::saveLevelSetVTK,
            "Write the level set to a ParaView VTK file.",
            py::arg("fileName"), py::arg("levelSet"), py::arg("isVelocity") = false,
            py::arg("isGradient") = false)

        .def("saveLevelSetTXT", (void (InputOutput::*)(const unsigned int&,
            const LevelSet&, const std::string&, bool) const) &InputOutput::saveLevelSetTXT,
            "Write the level set to a plain text file.",
            py::arg("datapoint"), py::arg("levelSet"), py::arg("outputDirectory") = "",
            py::arg("isXY") = false)

        .def("saveLevelSetBIN", (void (InputOutput::*)(const unsigned int&,
            const LevelSet&, const std::string&) const) &InputOutput::saveLevelSetBIN,
            "Write the level set to a binary file.",
            py::arg("datapoint"), py::arg("levelSet"), py::arg("outputDirectory") = "")

        .def("saveLevelSetTXT", (void (InputOutput::*)(const std::string&,
            const LevelSet&, bool) const) &InputOutput::saveLevelSetTXT,
            "Write the level set to a plain text file.",
            py::arg("filename"), py::arg("levelSet"), py::arg("isXY") = false)

        .def("saveLevelSetBIN", (void (InputOutput::*)(const std::string&,
            const LevelSet&) const) &InputOutput::saveLevelSetBIN,
            "Write the level set to a binary file.",
            py::arg("fileName"), py::arg("levelSet"))

        .def("loadLevelSetTXT", (void (InputOutput::*)(const unsigned int&,
            LevelSet&, const std::string&, bool) const) &InputOutput::loadLevelSetTXT,
            "Load the level set from a plain text file.",
            py::arg("datapoint"), py::arg("levelSet"), py::arg("outputDirectory") = "",
            py::arg("isXY") = false)

        .def("loadLevelSetTXT", (void (InputOutput::*)(const std::string&,
            LevelSet&, bool) const) &InputOutput::loadLevelSetTXT,
            "Load the level set from a plain text file.",
            py::arg("filename"), py::arg("levelSet"), py::arg("isXY") = false)

        .def("loadLevelSetBIN", (void (InputOutput::*)(const unsigned int&,
            LevelSet&, const std::string&) const) &InputOutput::loadLevelSetBIN,
            "Load the level from a binary file.",
            py::arg("datapoint"), py::arg("levelSet"), py::arg("outputDirectory") = "")

        .def("loadLevelSetBIN", (void (InputOutput::*)(const std::string&,
            LevelSet&) const) &InputOutput::loadLevelSetBIN,
            "Load the level set from a binary file.",
            py::arg("fileName"), py::arg("levelSet"))

        .def("saveBoundaryPointsTXT", (void (InputOutput::*)(const unsigned int&,
            const Boundary&, const std::string&) const) &InputOutput::saveBoundaryPointsTXT,
            "Save boundary point information to a plain text file.",
            py::arg("datapoint"), py::arg("boundary"), py::arg("outputDirectory") = "")

        .def("saveBoundaryPointsTXT", (void (InputOutput::*)(const std::string&,
            const Boundary&) const) &InputOutput::saveBoundaryPointsTXT,
            "Save boundary point information to a plain text file.",
            py::arg("fileName"), py::arg("boundary"))

        .def("saveBoundarySegmentsTXT", (void (InputOutput::*)(const unsigned int&,
            const Boundary&, const std::string&) const) &InputOutput::saveBoundarySegmentsTXT,
            "Save boundary segment information to a plain text file.",
            py::arg("datapoint"), py::arg("boundary"), py::arg("outputDirectory") = "")

        .def("saveBoundarySegmentsTXT", (void (InputOutput::*)(const std::string&,
            const Boundary&) const) &InputOutput::saveBoundarySegmentsTXT,
            "Save boundary segment information to a plain text file.",
            py::arg("fileName"), py::arg("boundary"))

        .def("saveAreaFractionsVTK", (void (InputOutput::*)(const unsigned int&,
            const Mesh&, const std::string&) const) &InputOutput::saveAreaFractionsVTK,
            "Write the element area fractions to a ParaView VTK file.",
            py::arg("datapoint"), py::arg("mesh"), py::arg("outputDirectory") = "")

        .def("saveAreaFractionsVTK", (void (InputOutput::*)(const std::string&,
            const Mesh&) const) &InputOutput::saveAreaFractionsVTK,
            "Write the element area fractions to a ParaView VTK file.",
            py::arg("fileName"), py::arg("mesh"))

        .def("saveAreaFractionsTXT", (void (InputOutput::*)(const unsigned int&,
            const Mesh&, const std::string&, bool) const) &InputOutput::saveAreaFractionsTXT,
            "Write the element area fractions to a plain text file.",
            py::arg("datapoint"), py::arg("mesh"), py::arg("outputDirectory") = "",
            py::arg("isXY") = false)

        .def("saveAreaFractionsTXT", (void (InputOutput::*)(const std::string&,
            const Mesh&, bool) const) &InputOutput::saveAreaFractionsTXT,
            "Write the element area fractions to a plain text file.",
            py::arg("fileName"), py::arg("mesh"), py::arg("isXY") = false);
}
