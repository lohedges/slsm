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
        // Start and end points of the segment.
        unsigned int start = boundary.segments[i].start;
        unsigned int end = boundary.segments[i].end;

        // Coordinates.
        double x, y;

        // First point.

        // Start is a boundary point.
        if (start >= mesh.nNodes)
        {
            start -= mesh.nNodes;
            x = boundary.points[start].x;
            y = boundary.points[start].y;
        }
        else
        {
            x = mesh.nodes[start].coord.x;
            y = mesh.nodes[start].coord.y;
        }

        // Write boundary point to file.
        fprintf(pFile, "%lf %lf\n", x, y);

        // Second point.

        // End is a boundary point.
        if (end >= mesh.nNodes)
        {
            end -= mesh.nNodes;
            x = boundary.points[end].x;
            y = boundary.points[end].y;
        }
        else
        {
            x = mesh.nodes[end].coord.x;
            y = mesh.nodes[end].coord.y;
        }

        // Write boundary point to file.
        fprintf(pFile, "%lf %lf\n\n", x, y);
    }

    fclose(pFile);

    return;

error:
    exit(EXIT_FAILURE);
}
