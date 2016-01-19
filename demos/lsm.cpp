/*
  Copyright (c) 2015-2016 Lester Hedges <lester.hedges+lsm@gmail.com>

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

#include <iostream>
#include <fstream>

#include "lsm.h"

int main(int argc, char** argv)
{
    // Print git commit info, if present.
#ifdef COMMIT
    std::cout << "Git commit: " << COMMIT << "\n";
#endif

    // Print git branch info, if present.
#ifdef COMMIT
    std::cout << "Git branch: " << BRANCH << "\n";
#endif

    // Initialise a mesh.
    //Mesh mesh(160, 80, false);
    Mesh mesh(100, 100, false);

    // Create a hole.
    std::vector<Hole> holes;
    holes.push_back(Hole(50, 50, 20));

    // Initialise the level set function (from hole vector).
    LevelSet levelSet(mesh, 3, holes);

    // Initialise the level set function (default Swiss cheese).
    //LevelSet levelSet(mesh, 3);

    // Initialise io object.
    InputOutput io;

    // Read Peter's level set.
    /*std::ifstream infile("input.txt");
    unsigned int i = 0;
    while (infile >> levelSet.signedDistance[i])
        i++;*/

    // Re-initialise the level set to a signed distance function.
    levelSet.reinitialise();

    // Initialise the boundary object.
    Boundary boundary(mesh, levelSet);

    // Discretise the boundary.
    boundary.discretise();

    // Compute element areas.
    boundary.computeAreaFractions();

    // Compute number of holes.
    unsigned int nHoles = boundary.computeHoles();

    // Print some statistics.
    std::cout << "\nBoundary length:   " << boundary.length << '\n';
    std::cout << "Material fraction: " << (boundary.area / (mesh.width*mesh.height)) << '\n';
    std::cout << "Number of holes:   " << nHoles << '\n';

    // Save boundary points and segments (txt file).
    io.saveBoundaryPointsTXT(1, boundary);
    io.saveBoundarySegmentsTXT(1, mesh, boundary);

    // Save element area fractions.
    io.saveAreaFractionsVTK(1, mesh);
    io.saveAreaFractionsTXT(1, mesh, "", true);

    for (unsigned int i=0;i<boundary.points.size();i++)
        boundary.points[i].sensitivities[0] = -1.0;

    // Data structures for Optimise class.
    std::vector<double> constraintDistances, lambdas;
    double timeStep;

    // Resize vectors.
    constraintDistances.resize(1);
    lambdas.resize(1);

    // Test optimisation class.
    Optimise optimise(boundary.points, constraintDistances, lambdas, timeStep);
    double areaChange = optimise.solve();

    // Print optimisation results.
    std::cout << "Time step:         " << timeStep << '\n';
    std::cout << "Area change:       " << areaChange << '\n';

    // Extend boundary point velocities to all level set nodes.
    levelSet.computeVelocities(boundary.points);

    // Save LSF info (ParaView and txt file).
    io.saveLevelSetVTK(1, mesh, levelSet);
    io.saveLevelSetTXT(1, mesh, levelSet, "", true);

    return (EXIT_SUCCESS);
}
