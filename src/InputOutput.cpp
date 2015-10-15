/*
 * File:	InputOutput.cpp
 * Author:	lester
 */

#include "InputOutput.h"

InputOutput::InputOutput() {}

void InputOutput::saveLevelSetVTK(const unsigned int& datapoint,
    const Mesh& mesh, const LevelSet& levelSet, const std::string& outputDirectory) const
{
    std::ostringstream fileName, num;

    num.str("");
    num.width(4);
    num.fill('0');
    num << std::right << datapoint;

    fileName.str("");
    if (!outputDirectory.empty()) fileName << outputDirectory << "/";
    fileName << "lsf_" << num.str() << ".vtk";

    saveLevelSetVTK(fileName, mesh, levelSet);
}

void InputOutput::saveLevelSetVTK(const std::ostringstream& fileName, const Mesh& mesh, const LevelSet& levelSet) const
{
    FILE *pFile;

    pFile = fopen(fileName.str().c_str(), "w");

    // Set up paraview header information.
    fprintf(pFile, "# vtk DataFile Version 3.0\n");
    fprintf(pFile, "Para0\n");
    fprintf(pFile, "ASCII\n");
    fprintf(pFile, "DATASET RECTILINEAR_GRID\n");
    fprintf(pFile, "DIMENSIONS %d %d %d\n", 1 + mesh.width, 1 + mesh.height, 1);
    fprintf(pFile, "X_COORDINATES %d int\n", 1 + mesh.width);
    for (unsigned int i=0;i<=mesh.width;i++) {fprintf(pFile, "%d ", i);}
    fprintf(pFile, "\nY_COORDINATES %d int\n", 1 + mesh.height);
    for (unsigned int i=0;i<=mesh.height;i++) {fprintf(pFile, "%d ", i);}
    fprintf(pFile, "\nZ_COORDINATES 1 int\n0\n\n");

    // Print out the nodal signed distance information.
    fprintf(pFile, "POINT_DATA %d\n", mesh.nNodes);
    fprintf(pFile, "SCALARS level-set float 1\n");
    fprintf(pFile, "LOOKUP_TABLE default\n");
    for (unsigned int i=0;i<mesh.nNodes;i++)
        fprintf(pFile, "%lf\n", levelSet.signedDistance[i]);

    fclose(pFile);
}
