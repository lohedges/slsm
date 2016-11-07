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

// Make sure we have PI defined.
#ifndef M_PI
    #define M_PI 3.1415926535897932384626433832795
#endif

#include <fstream>
#include <iostream>

#include "slsm.h"

/*! \file minimise_perimeter.cpp

    \brief An example code showing a hole shrinking during unconstrained
    perimeter minimisation.

    This demo provides a simple test to confirm that a single hole will
    shrink with velocity proportional to its curvature.

    At each iteration we optimise for the velocity vector that maximises
    the reduction in the boundary perimeter, subject to the CFL displacement
    limit. While this trivial optimisation could be performed by hand, it
    provides a useful test of our numerical implementation, especially the
    scaling that is required for stability inside the Optimise class.

    For a perfect continuous circle, all boundary points should move at a
    velocity proportional to the curvature, 1 / R. The time step will be
    rescaled if the curvature is too large and boundary point velocities
    cause violation of the CFL condition.

    We compute the velocity by measuring the distance that the hole boundary
    displaces as a function of time. This can be done by measuring the change
    in the boundary length over time and comparing the radius at each time
    interval (radius = perimeter / (2 x pi)).

    The output file, "minimise_perimeter.txt", contains the measured distance vs
    time data for the optmisation run. Additional columns provide data for the
    velocity and mean curvature at each time step (computed by the code) as well
    as the analytical curvature estimate of 1 / R. Level set information for each
    sample interval is written to ParaView readable VTK files, "level-set_*.vtk".
    Boundary segment data is written to "boundary-segments_*.txt".
 */

int main(int argc, char** argv)
{
    // Print git commit info, if present.
#ifdef COMMIT
    std::cout << "Git commit: " << COMMIT << "\n";
#endif

    // Print git branch info, if present.
#ifdef BRANCH
    std::cout << "Git branch: " << BRANCH << "\n";
#endif

    // Maximum displacement per iteration, in units of the mesh spacing.
    // This is the CFL limit.
    double moveLimit = 0.05;

    // Set maximum running time.
    double maxTime = 3000;

    // Set sampling interval.
    double sampleInterval = 30;

    // Set time of the next sample.
    double nextSample = 30;

    // Create a hole at position (100, 100) with a radius of 80 grid units.
    std::vector<slsm::Hole> holes;
    holes.push_back(slsm::Hole(100, 100, 80));

    // Initialise a 200x200 level set domain.
    slsm::LevelSet levelSet(200, 200, holes, moveLimit, 6, true);

    // Initialise io object.
    slsm::InputOutput io;

    // Reinitialise the level set to a signed distance function.
    levelSet.reinitialise();

    // Initialise the boundary object.
    slsm::Boundary boundary;

    // Perform initial boundary discretisation.
    boundary.discretise(levelSet);

    // Compute the initial element area fractions.
    levelSet.computeAreaFractions(boundary);

    // Compute the initial boundary point normal vectors.
    boundary.computeNormalVectors(levelSet);

    // Number of cycles since signed distance reinitialisation.
    unsigned int nReinit = 0;

    // Running time.
    double time = 0;

    // Time measurements.
    std::vector<double> times;

    // Boundary length measurements.
    std::vector<double> lengths;

    // Boundary curvature measurements.
    std::vector<double> curvatures;

    /* Lambda values for the optimiser.
       These are reused, i.e. the solution from the current iteration is
       used as an estimate for the next, hence we declare the vector
       outside of the main loop.
     */
    std::vector<double> lambdas(1);

    std::cout << "\nStarting unconstrained perimeter minimisation demo...\n\n";

    // Print output header.
    printf("--------------------------\n");
    printf("%6s %8s %10s\n", "Time", "Length", "Curvature");
    printf("--------------------------\n");

    // Integrate until we exceed the maximum time.
    while (time < maxTime)
    {
        // Zero the mean curvature.
        double curvature = 0;

        // Initialise the sensitivity object.
        slsm::Sensitivity sensitivity;

        // Initialise the sensitivity callback.
        using namespace std::placeholders;
        slsm::SensitivityCallback callback = std::bind(&slsm::Boundary::computePerimeter, boundary, _1);

        // Assign boundary point sensitivities.
        for (unsigned int i=0;i<boundary.points.size();i++)
        {
            boundary.points[i].sensitivities[0] =
                sensitivity.computeSensitivity(boundary.points[i], callback);

            curvature += boundary.points[i].sensitivities[0];
        }

        // Compute mean curvature.
        curvature /= boundary.points.size();

        // Time step associated with the iteration.
        double timeStep;

        /* Initialise the optimisation object.

           The Optimise class is a lightweight object so there is no cost for
           reinitialising at every iteration. A smart compiler will optimise
           this anyway, i.e. the same memory space will be reused. It is better
           to place objects in the correct scope in order to aid readability
           and to avoid unintended name clashes, etc.

           Since there are no constraints we pass an empty vector for the
           constraint distances argument.
         */
        slsm::Optimise optimise(boundary.points, std::vector<double>(),
            lambdas, timeStep, levelSet.moveLimit, false);

        // Perform the optimisation.
        optimise.solve();

        // Extend boundary point velocities to all narrow band nodes.
        levelSet.computeVelocities(boundary.points);

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
        boundary.discretise(levelSet);

        // Compute the element area fractions.
        levelSet.computeAreaFractions(boundary);

        // Compute the boundary point normal vectors.
        boundary.computeNormalVectors(levelSet);

        // Increment the time.
        time += timeStep;

        // Check if the next sample time has been reached.
        while (time >= nextSample)
        {
            // Record the time, boundary length, and mean curvature.
            times.push_back(time);
            lengths.push_back(boundary.length);
            curvatures.push_back(curvature);

            // Update the time of the next sample.
            nextSample += sampleInterval;

            // Print statistics.
            printf("%6.1f %8.1f %10.4f\n", time, boundary.length, curvature);

            // Write level set and boundary segments to file.
            io.saveLevelSetVTK(times.size(), levelSet);
            io.saveBoundarySegmentsTXT(times.size(), boundary);
        }
    }

    // Distance measurements.
    std::vector<double> distances;

    // Compute the distance moved at each time interval.
    for (unsigned int i=0;i<times.size();i++)
        distances.push_back((lengths[0] - lengths[i]) / (2 * M_PI));

    // Print results to file.
    FILE *pFile;
    pFile = fopen("minimise_perimeter.txt", "w");
    for (unsigned int i=1;i<times.size();i++)
    {
        // Distance and time increments.
        double deltaDist = distances[i] - distances[i-1];
        double deltaTime = times[i] - times[i-1];

        fprintf(pFile, "%lf %lf %lf %lf %lf\n", times[i] - times[0],
            distances[i], deltaDist / deltaTime, curvatures[i], ((2 * M_PI) / lengths[i]));
    }
    fclose(pFile);

    std::cout << "\nDone!\n";

    return (EXIT_SUCCESS);
}
