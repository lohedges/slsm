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
    lsm::Mesh mesh(200, 200, false);

    // Create a hole.
    std::vector<lsm::Hole> holes;
    holes.push_back(lsm::Hole(100, 100, 80));

    // Initialise the level set function (from hole vector).
    lsm::LevelSet levelSet(mesh, holes, 0.5, 3);

    // Initialise io object.
    lsm::InputOutput io;

    // Reinitialise the level set to a signed distance function.
    levelSet.reinitialise();

    // Initialise the boundary object.
    lsm::Boundary boundary(mesh, levelSet);

    // Number of cycles since reinitialisation.
    unsigned int nReinit = 0;

    // Integrate for 100 time steps.
    for (unsigned int i=0;i<100;i++)
    {
        std::cout << "\nStarting iteration: " << i+1 << '\n';

        // Discretise the boundary.
        boundary.discretise();

        // Compute element areas.
        boundary.computeAreaFractions();

        // Print some statistics.
        std::cout << "Boundary length:    " << boundary.length << '\n';
        std::cout << "Material fraction:  " << (boundary.area / (mesh.width*mesh.height)) << '\n';

        // Fill sensitivities with dummy values.
        for (unsigned int i=0;i<boundary.points.size();i++)
            boundary.points[i].sensitivities[0] = 1.0;

        // Data structures for Optimise class.
        std::vector<double> constraintDistances, lambdas;
        double timeStep;

        // Resize vectors.
        constraintDistances.resize(1);
        lambdas.resize(1);

        // Test optimisation class.
        lsm::Optimise optimise(boundary.points, constraintDistances, lambdas, timeStep);
        double areaChange = optimise.solve();

        // Print optimisation results.
        std::cout << "Time step:          " << timeStep << '\n';
        std::cout << "Area change:        " << areaChange << '\n';

        // Extend boundary point velocities to all level set nodes.
        levelSet.computeVelocities(boundary.points);

        // Compute gradient of signed distance function.
        levelSet.computeGradients(timeStep);

        // Save LSF info (ParaView and txt file).
        io.saveLevelSetVTK(i+1, mesh, levelSet);
        io.saveLevelSetTXT(i+1, mesh, levelSet, "", true);

        // Save boundary points and segments (txt file).
        io.saveBoundaryPointsTXT(i+1, boundary);
        io.saveBoundarySegmentsTXT(i+1, mesh, boundary);

        // Save element area fractions.
        io.saveAreaFractionsVTK(i+1, mesh);
        io.saveAreaFractionsTXT(i+1, mesh, "", true);

        // Update the level set function.
        bool isReinitialised = levelSet.update(timeStep);

        // Reinitialise the signed distance function.
        if (!isReinitialised)
        {
            // Reinitialise every 20 iterations.
            if (nReinit == 20)
            {
                levelSet.reinitialise();
                nReinit = 0;
            }
        }
        else nReinit = 0;

        // Increment number of steps since reinitialisation.
        nReinit++;
    }

    return (EXIT_SUCCESS);
}
