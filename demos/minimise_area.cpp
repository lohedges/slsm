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

#include "lsm.h"

/*! \file minimise_area.cpp

    \brief An example code showing a hole shrinking during unconstrained
    area maximisation.

    This demo provides a simple test to confirm that a single hole will
    shrink with constant velocity when we perform unconstrained optmisation,
    i.e. maximise the material area fraction subject to no constraint.

    At each iteration we optimise for the velocity vector that maximises
    the change in the material area fraction, subject to the CFL displacement
    limit. While this trivial optimisation could be performed by hand, it
    provides a useful test of our numerical implementation, especially the
    scaling that is required for stability inside the Optimise class.

    All boundary points should move at a constant unit velocity, where the time
    step is dictated by the CFL limit. For example, if the CFL limit is 0.5
    (half a grid spacing) then boundary points moving at unit velocity will
    have an associated time step of 0.5, i.e velocity x time = CFL.

    We test that the hole shrinks with unit velocity by measuring the distance
    that the hole boundary displaces as a function of time. This can be done
    by measuring the change in the boundary length over time and comparing the
    radius at each time interval (radius = perimeter / (2 x pi)).

    The output file, "minimise_area.txt", contains the measured distance vs
    time data for the optmisation run. Level set information for each
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
    double moveLimit = 0.5;

    // Set maximum running time.
    double maxTime = 50;

    // Set sampling interval.
    double sampleInterval = 1;

    // Set time of the next sample.
    double nextSample = 1;

    // Initialise a 200x200 non-periodic mesh.
    lsm::Mesh mesh(200, 200, false);

    // Create a hole in the centre with a radius of 80 grid units.
    std::vector<lsm::Hole> holes;
    holes.push_back(lsm::Hole(100, 100, 80));

    // Initialise the level set object (from the hole vector).
    lsm::LevelSet levelSet(mesh, holes, moveLimit, 6, true);

    // Initialise io object.
    lsm::InputOutput io;

    // Reinitialise the level set to a signed distance function.
    levelSet.reinitialise();

    // Initialise the boundary object.
    lsm::Boundary boundary(levelSet);

    // Perform initial boundary discretisation.
    boundary.discretise();

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
    std::vector<double> lambdas(1);

    std::cout << "\nStarting unconstrained area minimisation demo...\n\n";

    // Print output header.
    printf("---------------\n");
    printf("%6s %8s\n", "Time", "Length");
    printf("---------------\n");

    // Integrate until we exceed the maximum time.
    while (time < maxTime)
    {
        // Assign boundary point sensitivities.
        for (unsigned int i=0;i<boundary.points.size();i++)
            boundary.points[i].sensitivities[0] = 1.0;

        // Time step associated with the iteration.
        double timeStep;

        // Constraint distance vector. Since there are no constraints the
        // vector can be left unassigned.
        std::vector<double> constraintDistances(1);

        /* Initialise the optimisation object.

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
        boundary.discretise();

        // Increment the time.
        time += timeStep;

        // Check if the next sample time has been reached.
        while (time >= nextSample)
        {
            // Record the time and boundary length.
            times.push_back(time);
            lengths.push_back(boundary.length);

            // Update the time of the next sample.
            nextSample += sampleInterval;

            // Print statistics.
            printf("%6.1f %8.1f\n", time, boundary.length);

            // Write level set and boundary segments to file.
            io.saveLevelSetVTK(times.size(), levelSet);
            io.saveBoundarySegmentsTXT(times.size(), boundary);
        }
    }

    // Print results to file (distance vs time).
    FILE *pFile;
    pFile = fopen("minimise_area.txt", "w");
    for (unsigned int i=0;i<times.size();i++)
        fprintf(pFile, "%lf %lf\n", times[i] - times[0], (lengths[0] - lengths[i]) / (2 * M_PI));
    fclose(pFile);

    std::cout << "\nDone!\n";

    return (EXIT_SUCCESS);
}
