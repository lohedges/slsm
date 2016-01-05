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
    Mesh mesh(160, 80, false);

    // Create a hole.
    /*std::vector<Hole> holes;
    holes.push_back(Hole(5, 5, 5));

    // Initialise the level set function (from hole vector).
    LevelSet levelSet(mesh, 3, holes);*/

    // Initialise the level set function (default Swiss cheese).
    LevelSet levelSet(mesh, 3);

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
    std::cout << "Boundary length:   " << boundary.length << '\n';
    std::cout << "Material fraction: " << (boundary.area / (mesh.width*mesh.height)) << '\n';
    std::cout << "Number of holes:   " << nHoles << '\n';

    // Save LSF info (ParaView and txt file).
    io.saveLevelSetVTK(1, mesh, levelSet);
    io.saveLevelSetTXT(1, mesh, levelSet, "", true);

    // Save boundary points and segments (txt file).
    io.saveBoundaryPointsTXT(1, boundary);
    io.saveBoundarySegmentsTXT(1, mesh, boundary);

    // Save element area fractions.
    io.saveAreaFractionsVTK(1, mesh);
    io.saveAreaFractionsTXT(1, mesh, "", true);

    // Test optimisation class.
    std::vector<double> tmp1, tmp2, tmp3;
    tmp1.resize(boundary.points.size());
    tmp2.resize(1);
    tmp3.resize(2);
    Optimise optimise(boundary.points, tmp2, tmp3, tmp1);

    return (EXIT_SUCCESS);
}
