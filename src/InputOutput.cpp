/*
  Copyright (c) 2015-2016 Lester Hedges <lester.hedges+slsm@gmail.com>

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

#include "Boundary.h"
#include "Debug.h"
#include "InputOutput.h"
#include "LevelSet.h"
#include "Mesh.h"

/*! \file InputOutput.cpp
    \brief A class for reading and writing data.
 */

namespace slsm
{
    InputOutput::InputOutput() {}

    void InputOutput::saveLevelSetVTK(const unsigned int& datapoint, const LevelSet& levelSet,
        bool isVelocity, bool isGradient, const std::string& outputDirectory) const
    {
        std::ostringstream fileName, num;

        num.str("");
        num.width(4);
        num.fill('0');
        num << std::right << datapoint;

        fileName.str("");
        if (!outputDirectory.empty()) fileName << outputDirectory << "/";
        fileName << "level-set_" << num.str() << ".vtk";

        saveLevelSetVTK(fileName, levelSet);
    }

    void InputOutput::saveLevelSetVTK(const std::ostringstream& fileName,
        const LevelSet& levelSet, bool isVelocity, bool isGradient) const
    {
        FILE *pFile;

        pFile = fopen(fileName.str().c_str(), "w");

        if (pFile == NULL)
            slsm_sentinel("Write error, cannot open file %s", fileName.str().c_str())

        // Set up ParaView header information.
        fprintf(pFile, "# vtk DataFile Version 3.0\n");
        fprintf(pFile, "Para0\n");
        fprintf(pFile, "ASCII\n");
        fprintf(pFile, "DATASET RECTILINEAR_GRID\n");
        fprintf(pFile, "DIMENSIONS %d %d %d\n", 1 + levelSet.mesh.width, 1 + levelSet.mesh.height, 1);
        fprintf(pFile, "X_COORDINATES %d int\n", 1 + levelSet.mesh.width);
        for (unsigned int i=0;i<=levelSet.mesh.width;i++) {fprintf(pFile, "%d ", i);}
        fprintf(pFile, "\nY_COORDINATES %d int\n", 1 + levelSet.mesh.height);
        for (unsigned int i=0;i<=levelSet.mesh.height;i++) {fprintf(pFile, "%d ", i);}
        fprintf(pFile, "\nZ_COORDINATES 1 int\n0\n\n");
        fprintf(pFile, "POINT_DATA %d\n", levelSet.mesh.nNodes);

        // Write the nodal signed distance to file.
        fprintf(pFile, "SCALARS distance float 1\n");
        fprintf(pFile, "LOOKUP_TABLE default\n");
        for (unsigned int i=0;i<levelSet.mesh.nNodes;i++)
            fprintf(pFile, "%lf\n", levelSet.signedDistance[i]);

        // Write the nodal velocity to file.
        if (isVelocity)
        {
            fprintf(pFile, "SCALARS velocity float 1\n");
            fprintf(pFile, "LOOKUP_TABLE default\n");
            for (unsigned int i=0;i<levelSet.mesh.nNodes;i++)
                fprintf(pFile, "%lf\n", levelSet.velocity[i]);
        }

        // Write the nodal gradient to file.
        if (isGradient)
        {
            fprintf(pFile, "SCALARS gradient float 1\n");
            fprintf(pFile, "LOOKUP_TABLE default\n");
            for (unsigned int i=0;i<levelSet.mesh.nNodes;i++)
                fprintf(pFile, "%lf\n", levelSet.gradient[i]);
        }

        fclose(pFile);

        return;

    error:
        exit(EXIT_FAILURE);
    }

    void InputOutput::saveLevelSetTXT(const unsigned int& datapoint,
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

        saveLevelSetTXT(fileName, levelSet, isXY);
    }

    void InputOutput::saveLevelSetTXT(const std::ostringstream& fileName,
        const LevelSet& levelSet, bool isXY) const
    {
        FILE *pFile;

        pFile = fopen(fileName.str().c_str(), "w");

        if (pFile == NULL)
            slsm_sentinel("Write error, cannot open file %s", fileName.str().c_str())

        // Write the nodal signed distance to file.
        for (unsigned int i=0;i<levelSet.mesh.nNodes;i++)
        {
            if (isXY) fprintf(pFile, "%lf %lf ", levelSet.mesh.nodes[i].coord.x, levelSet.mesh.nodes[i].coord.y);
            fprintf(pFile, "%lf %lf %lf\n", levelSet.signedDistance[i], levelSet.velocity[i], levelSet.gradient[i]);
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
            slsm_sentinel("Write error, cannot open file %s", fileName.str().c_str())

        // Write the boundary points to file.
        for (unsigned int i=0;i<boundary.nPoints;i++)
            fprintf(pFile, "%lf %lf %lf\n",
                boundary.points[i].coord.x, boundary.points[i].coord.y, boundary.points[i].length);

        fclose(pFile);

        return;

    error:
        exit(EXIT_FAILURE);
    }

    void InputOutput::saveBoundarySegmentsTXT(const unsigned int& datapoint,
        const Boundary& boundary, const std::string& outputDirectory) const
    {
        std::ostringstream fileName, num;

        num.str("");
        num.width(4);
        num.fill('0');
        num << std::right << datapoint;

        fileName.str("");
        if (!outputDirectory.empty()) fileName << outputDirectory << "/";
        fileName << "boundary-segments_" << num.str() << ".txt";

        saveBoundarySegmentsTXT(fileName, boundary);
    }

    void InputOutput::saveBoundarySegmentsTXT(const std::ostringstream& fileName, const Boundary& boundary) const
    {
        FILE *pFile;

        pFile = fopen(fileName.str().c_str(), "w");

        if (pFile == NULL)
            slsm_sentinel("Write error, cannot open file %s", fileName.str().c_str())

        // Write the boundary points to file.
        for (unsigned int i=0;i<boundary.nSegments;i++)
        {
            // Start and end points of the segment.
            unsigned int start = boundary.segments[i].start;
            unsigned int end = boundary.segments[i].end;

            // Coordinates.
            double x, y;

            // First point.
            x = boundary.points[start].coord.x;
            y = boundary.points[start].coord.y;

            // Write boundary point to file.
            fprintf(pFile, "%lf %lf\n", x, y);

            // Second point.
            x = boundary.points[end].coord.x;
            y = boundary.points[end].coord.y;

            // Write boundary point to file.
            fprintf(pFile, "%lf %lf\n\n", x, y);
        }

        fclose(pFile);

        return;

    error:
        exit(EXIT_FAILURE);
    }

    void InputOutput::saveAreaFractionsVTK(const unsigned int& datapoint,
        const Mesh& mesh, const std::string& outputDirectory) const
    {
        std::ostringstream fileName, num;

        num.str("");
        num.width(4);
        num.fill('0');
        num << std::right << datapoint;

        fileName.str("");
        if (!outputDirectory.empty()) fileName << outputDirectory << "/";
        fileName << "area_" << num.str() << ".vtk";

        saveAreaFractionsVTK(fileName, mesh);
    }

    void InputOutput::saveAreaFractionsVTK(const std::ostringstream& fileName, const Mesh& mesh) const
    {
        FILE *pFile;

        pFile = fopen(fileName.str().c_str(), "w");

        if (pFile == NULL)
            slsm_sentinel("Write error, cannot open file %s", fileName.str().c_str())

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

        // Write the element area fractions to file.
        fprintf(pFile, "CELL_DATA %d\n", mesh.nElements);
        fprintf(pFile, "SCALARS area float 1\n");
        fprintf(pFile, "LOOKUP_TABLE default\n");
        for (unsigned int i=0;i<mesh.nElements;i++)
            fprintf(pFile, "%lf\n", mesh.elements[i].area);

        fclose(pFile);

        return;

    error:
        exit(EXIT_FAILURE);
    }

    void InputOutput::saveAreaFractionsTXT(const unsigned int& datapoint,
        const Mesh& mesh, const std::string& outputDirectory, bool isXY) const
    {
        std::ostringstream fileName, num;

        num.str("");
        num.width(4);
        num.fill('0');
        num << std::right << datapoint;

        fileName.str("");
        if (!outputDirectory.empty()) fileName << outputDirectory << "/";
        fileName << "area_" << num.str() << ".txt";

        saveAreaFractionsTXT(fileName, mesh, isXY);
    }

    void InputOutput::saveAreaFractionsTXT(const std::ostringstream& fileName, const Mesh& mesh, bool isXY) const
    {
        FILE *pFile;

        pFile = fopen(fileName.str().c_str(), "w");

        if (pFile == NULL)
            slsm_sentinel("Write error, cannot open file %s", fileName.str().c_str())

        // Write the element area fractions to file.
        for (unsigned int i=0;i<mesh.nElements;i++)
        {
            if (isXY) fprintf(pFile, "%lf %lf ", mesh.elements[i].coord.x, mesh.elements[i].coord.y);
            fprintf(pFile, "%lf\n", mesh.elements[i].area);
        }

        fclose(pFile);

        return;

    error:
        exit(EXIT_FAILURE);
    }
}
