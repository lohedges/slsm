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

#include <fstream>
#include <iostream>

#include "lsm.h"

/*! \file minimise_perimeter.cpp
    \brief An example code showing perimeter minimisation with noise.

    The output file, "perimeter.txt", contains the measured boundary length vs
    time data for the optmisation run. Level set information for each sample
    interval is written to ParaView readable VTK files, "level-set_*.vtk".
    Boundary segment data is written to "boundary-segments_*.txt".
 */

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

    // Maximum displacement per iteration, in units of the mesh spacing.
    // This is the CFL limit.
    double moveLimit = 0.5;

    // Temperature of the thermal bath.
    double temperature = 0.02;

    // Set maximum running time.
    double maxTime = 1000;

    // Set minimum material area.
    double minArea = 0.6;

    // Set sampling interval.
    double sampleInterval = 10;

    // Set time of the next sample.
    double nextSample = 10;

    // Initialise a 200x200 non-periodic mesh.
    lsm::Mesh mesh(200, 200, false);

    // Calculate the area of the mesh.
    double meshArea = mesh.width * mesh.height;

    // Initialise the level set object.
    lsm::LevelSet levelSet(mesh, moveLimit, 6, true);

    // Create a solid slab of material with a small square in the middle.
    for (unsigned int i=0;i<mesh.nNodes;i++)
    {
        unsigned int x = int(mesh.nodes[i].coord.x);
        unsigned int y = int(mesh.nodes[i].coord.y);

        // Cut out square hole.
        if (x >= 80 && x <= 120 && y >= 80 && y <= 120)
            levelSet.signedDistance[i] = -1;

        // Fill with material.
        else levelSet.signedDistance[i] = 1;
    }

    // Initialise io object.
    lsm::InputOutput io;

    // Reinitialise the level set to a signed distance function.
    levelSet.reinitialise();

    // Initialise the boundary object.
    lsm::Boundary boundary(mesh, levelSet);

    // Perform initial boundary discretisation.
    boundary.discretise();
    boundary.computeAreaFractions();

    // Initialise random number generator.
    lsm::MersenneTwister rng;

    // Number of cycles since signed distance reinitialisation.
    unsigned int nReinit = 0;

    // Running time.
    double time = 0;

    // Time measurements.
    std::vector<double> times;

    // Boundary length measurements.
    std::vector<double> lengths;

    /* Lambda values for the optimiser.
       These are reused, i.e. the solution from the current iteration is
       used as an estimate for the next, hence we declare the vector
       outside of the main loop.
     */
    std::vector<double> lambdas(2);

    std::cout << "\nStarting perimeter minimisation demo...\n\n";

    // Print output header.
    printf("----------------------\n");
    printf("%6s %8s %6s\n", "Time", "Length", "Area");
    printf("----------------------\n");

    // Integrate until we exceed the maximim time.
    while (time < maxTime)
    {
        // Assign boundary point sensitivities.
        for (unsigned int i=0;i<boundary.points.size();i++)
        {
            boundary.points[i].sensitivities[0] = boundary.points[i].curvature;
            boundary.points[i].sensitivities[1] = -1.0;
        }

        // Time step associated with the iteration.
        double timeStep;

        // Constraint distance vector.
        std::vector<double> constraintDistances;

        // Push current distance from constraint violation into vector.
        constraintDistances.push_back(meshArea*minArea - boundary.area);

        /* Initialise the optimisation object for perimeter minimisation.

           The Optimise class is a lightweight object so there is no cost for
           reinitialising at every iteration. A smart compiler will optimise
           this anyway, i.e. the same memory space will be reused. It is better
           to place objects in the correct scope in order to aid readability
           and to avoid unintended name clashes, etc.
         */
        lsm::Optimise optimise(boundary.points, constraintDistances,
            lambdas, timeStep, levelSet.moveLimit, false);

        // Perform the optimisation.
        optimise.solve();

        // Extend boundary point velocities to all narrow band nodes.
        levelSet.computeVelocities(boundary.points, timeStep, temperature, rng);

        // Compute gradient of the signed distance function within the narrow band.
        levelSet.computeGradients();

        // Update the level set function.
        bool isReinitialised = levelSet.update(timeStep);

        // Reinitialise the signed distance function, if necessary.
        if (!isReinitialised)
        {
            // Reinitialise at least every 20 iterations.
            if (nReinit == 20)
            {
                levelSet.reinitialise();
                nReinit = 0;
            }
        }
        else nReinit = 0;

        // Increment the number of steps since reinitialisation.
        nReinit++;

        // Compute the new discretised boundary.
        boundary.discretise();
        boundary.computeAreaFractions();

        // Increment the time.
        time += timeStep;

        // Check if the next sample time has been reached.
        while (time >= nextSample)
        {
            // Record the time.
            times.push_back(time);

            // Update the time of the next sample.
            nextSample += sampleInterval;

            // Print statistics.
            printf("%6.2f %8.1f %6.2f\n", time, boundary.length, boundary.area / meshArea);

            // Write level set and boundary segments to file.
            io.saveLevelSetVTK(times.size(), mesh, levelSet);
            io.saveBoundarySegmentsTXT(times.size(), mesh, boundary);

            // Store the average radius.
            lengths.push_back(boundary.length);
        }
    }

    // Print results to file (perimeter vs time).
    FILE *pFile;
    pFile = fopen("perimeter.txt", "w");
    for (unsigned int i=0;i<times.size();i++)
        fprintf(pFile, "%lf %lf\n", times[i] - times[0], lengths[i]);

    std::cout << "\nDone!\n";

    return (EXIT_SUCCESS);
}
