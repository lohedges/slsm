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
    fileName << "level-set_" << num.str() << ".vtk";

    saveLevelSetVTK(fileName, mesh, levelSet);
}

void InputOutput::saveLevelSetVTK(const std::ostringstream& fileName, const Mesh& mesh, const LevelSet& levelSet) const
{
    FILE *pFile;

    pFile = fopen(fileName.str().c_str(), "w");

    if (pFile == NULL)
        sentinel("Write error, cannot open file %s", fileName.str().c_str())

    // Set up ParaView header information.
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

    // Write the nodal signed distance to file.
    fprintf(pFile, "POINT_DATA %d\n", mesh.nNodes);
    fprintf(pFile, "SCALARS level-set float 1\n");
    fprintf(pFile, "LOOKUP_TABLE default\n");
    for (unsigned int i=0;i<mesh.nNodes;i++)
        fprintf(pFile, "%lf\n", levelSet.signedDistance[i]);

    fclose(pFile);

    return;

error:
    exit(EXIT_FAILURE);
}

void InputOutput::saveLevelSetTXT(const unsigned int& datapoint, const Mesh& mesh,
    const LevelSet& levelSet, const std::string& outputDirectory, bool isXY) const
{
    std::ostringstream fileName, num;

    num.str("");
    num.width(4);
    num.fill('0');
    num << std::right << datapoint;

    fileName.str("");
    if (!outputDirectory.empty()) fileName << outputDirectory << "/";
    fileName << "level-set_" << num.str() << ".txt";

    saveLevelSetTXT(fileName, mesh, levelSet, isXY);
}

void InputOutput::saveLevelSetTXT(const std::ostringstream& fileName,
    const Mesh& mesh, const LevelSet& levelSet, bool isXY) const
{
    FILE *pFile;

    pFile = fopen(fileName.str().c_str(), "w");

    if (pFile == NULL)
        sentinel("Write error, cannot open file %s", fileName.str().c_str())

    // Write the nodal signed distance to file.
    for (unsigned int i=0;i<mesh.nNodes;i++)
    {
        if (isXY) fprintf(pFile, "%lf %lf ", mesh.nodes[i].coord.x, mesh.nodes[i].coord.y);
        fprintf(pFile, "%lf\n", levelSet.signedDistance[i]);
    }

    fclose(pFile);

    return;

error:
    exit(EXIT_FAILURE);
}

void InputOutput::saveBoundaryPointsTXT(const unsigned int& datapoint,
    const Boundary& boundary, const std::string& outputDirectory) const
{
    std::ostringstream fileName, num;

    num.str("");
    num.width(4);
    num.fill('0');
    num << std::right << datapoint;

    fileName.str("");
    if (!outputDirectory.empty()) fileName << outputDirectory << "/";
    fileName << "boundary-points_" << num.str() << ".txt";

    saveBoundaryPointsTXT(fileName, boundary);
}

void InputOutput::saveBoundaryPointsTXT(const std::ostringstream& fileName, const Boundary& boundary) const
{
    FILE *pFile;

    pFile = fopen(fileName.str().c_str(), "w");

    if (pFile == NULL)
        sentinel("Write error, cannot open file %s", fileName.str().c_str())

    // Write the boundary points to file.
    for (unsigned int i=0;i<boundary.nPoints;i++)
        fprintf(pFile, "%lf %lf\n", boundary.points[i].x, boundary.points[i].y);

    fclose(pFile);

    return;

error:
    exit(EXIT_FAILURE);
}

void InputOutput::saveBoundarySegmentsTXT(const unsigned int& datapoint,
    const Mesh& mesh, const Boundary& boundary, const std::string& outputDirectory) const
{
    std::ostringstream fileName, num;

    num.str("");
    num.width(4);
    num.fill('0');
    num << std::right << datapoint;

    fileName.str("");
    if (!outputDirectory.empty()) fileName << outputDirectory << "/";
    fileName << "boundary-segments_" << num.str() << ".txt";

    saveBoundarySegmentsTXT(fileName, mesh, boundary);
}

void InputOutput::saveBoundarySegmentsTXT(const std::ostringstream& fileName,
    const Mesh& mesh, const Boundary& boundary) const
{
    FILE *pFile;

    pFile = fopen(fileName.str().c_str(), "w");

    if (pFile == NULL)
        sentinel("Write error, cannot open file %s", fileName.str().c_str())

    // Write the boundary points to file.
    for (unsigned int i=0;i<boundary.nSegments;i++)
    {
        // Start and end nodes of the segment.
        unsigned int n1 = boundary.segments[i].node1;
        unsigned int n2 = boundary.segments[i].node2;

        // Coordinates.
        double x, y;

        // First point.

        // Node is a boundary point.
        if (n1 >= mesh.nNodes)
        {
            n1 -= mesh.nNodes;
            x = boundary.points[n1].x;
            y = boundary.points[n1].y;
        }
        else
        {
            x = mesh.nodes[n1].coord.x;
            y = mesh.nodes[n1].coord.y;
        }

        // Write boundary point to file.
        fprintf(pFile, "%lf %lf\n", x, y);

        // Second point.

        // Node is a boundary point.
        if (n2 >= mesh.nNodes)
        {
            n2 -= mesh.nNodes;
            x = boundary.points[n2].x;
            y = boundary.points[n2].y;
        }
        else
        {
            x = mesh.nodes[n2].coord.x;
            y = mesh.nodes[n2].coord.y;
        }

        // Write boundary point to file.
        fprintf(pFile, "%lf %lf\n", x, y);
    }

    fclose(pFile);

    return;

error:
    exit(EXIT_FAILURE);
}
