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
#include <iomanip>
#include <iostream>

#include "lsm.h"

/*! \file dumbbell.cpp

    \brief An example showing perimeter minimisation with a shape matching
    constraint.

    Here we construct a simple example system with two minima separated by
    a free energy barrier. The matched shape is a narrow-necked dumbbell
    constructed from two vertically offset, overlapping circles. The initial
    configuration is a circle centred in the upper lobe of the dumbbell.
    We construct two minima by modifying the perimeter (objective)
    sensitivities in the upper and lower half of the domain. Sensitivities
    are reduced in the lower half, hence it's possible to form a circle with
    a smaller perimeter at the same cost.

    To reach the global minimum in the lower lobe the shape must pass through
    the neck of the dumbbell. This transition requires a significant deformation
    to the interface and a corresponding increase in the perimeter of the zero
    contour. The pathway is not possible at zero temperature since it requires
    a fluctuation that is uphill in free energy. As such the circle remains
    trapped in the upper lobe.

    The output file, "dumbbell_*.txt", contains the measured
    perimeter and mismatch vs time data for the optmisation run. Level set
    information for each sample interval is written to ParaView readable VTK
    files, "level-set_*.vtk".  Boundary segment data is written to
    "boundary-segments_*.txt".
 */

// Sensitivity function prototype.
double computeSensitivity(const lsm::Coord&, const lsm::Mesh&, const lsm::LevelSet&);

// Mismatch function prototype.
double computeMismatch(const lsm::Mesh&, const std::vector<double>&);

// Perimeter function prototype.
double computePerimeter(const std::vector<lsm::BoundaryPoint>&, double, double);

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

    // Default temperature of the thermal bath.
    double temperature = 0.15;

    // Override temperature if command-line argument is passed.
    if (argc == 2) temperature = atof(argv[1]);

    // Set maximum running time.
    double maxTime = 6000;

    // Set maximumum area mismatch.
    double maxMismatch = 0.2;

    // Set sensitivity reduction factor.
    double reduce = 0.67;

    // Set sampling interval.
    double sampleInterval = 30;

    // Set time of the next sample.
    double nextSample = 30;

    // Initialise a 100x100 non-periodic mesh.
    lsm::Mesh mesh(100, 100, false);

    // Store the mesh area.
    double meshArea = mesh.width * mesh.height;

    // Half mesh height.
    double halfHeight = 0.5 * mesh.height;

    std::vector<lsm::Hole> initialHoles;
    std::vector<lsm::Hole> targetHoles;

    // Create a dumbbell from two vertically offset holes.
    targetHoles.push_back(lsm::Hole(50, 31, 20));
    targetHoles.push_back(lsm::Hole(50, 69, 20));

    // Initialise the system with a circle in the upper lobe.
    initialHoles.push_back(lsm::Hole(50, 69, 20));

    // Initialise the level set object.
    lsm::LevelSet levelSet(mesh, initialHoles, targetHoles, moveLimit, 6, true);

    // Initialise io object.
    lsm::InputOutput io;

    // Reinitialise the level set to a signed distance function.
    levelSet.reinitialise();

    // Initialise the boundary object.
    lsm::Boundary boundary(levelSet);

    // Initialise target area fraction vector.
    std::vector<double> targetArea(mesh.nElements);

    // Discretise the target structure.
    boundary.discretise(true);
    io.saveBoundarySegmentsTXT(0, boundary);

    // Compute the element area fractions.
    boundary.computeAreaFractions();

    // Store the target area fractions.
    for (unsigned int i=0;i<mesh.nElements;i++)
        targetArea[i] = mesh.elements[i].area;

    // Perform initial boundary discretisation.
    boundary.discretise();

    // Compute the element area fractions.
    boundary.computeAreaFractions();

    // Compute the initial boundary point normal vectors.
    boundary.computeNormalVectors();

    // Initialise random number generator.
    lsm::MersenneTwister rng;

    // Number of cycles since signed distance reinitialisation.
    unsigned int nReinit = 0;

    // Running time.
    double time = 0;

    // Time measurements.
    std::vector<double> times;

    // Boundary length measurements (objective).
    std::vector<double> lengths;

    // Area mismatch measurements (constraint).
    std::vector<double> mismatches;

    /* Lambda values for the optimiser.
       These are reused, i.e. the solution from the current iteration is
       used as an estimate for the next, hence we declare the vector
       outside of the main loop.
     */
    std::vector<double> lambdas(2);

    std::cout << "\nStarting dumbbell demo...\n\n";

    // Print output header.
    printf("--------------------------\n");
    printf("%6s %8s %10s\n", "Time", "Length", "Mismatch");
    printf("--------------------------\n");

    // Integrate until we exceed the maximum time.
    while (time < maxTime)
    {
        // Initialise the sensitivity object.
        lsm::Sensitivity sensitivity;

        // Initialise the sensitivity callback.
        using namespace std::placeholders;
        lsm::SensitivityCallback callback = std::bind(&lsm::Boundary::computePerimeter, boundary, _1);

        // Assign boundary point sensitivities.
        for (unsigned int i=0;i<boundary.points.size();i++)
        {
            boundary.points[i].sensitivities[0] =
                sensitivity.computeSensitivity(boundary.points[i], callback);
            boundary.points[i].sensitivities[1] =
                computeSensitivity(boundary.points[i].coord, mesh, levelSet);

            // Reduce perimeter sensitivities in the lower half.
            if (boundary.points[i].coord.y < halfHeight)
                boundary.points[i].sensitivities[0] *= reduce;
        }

        // Time step associated with the iteration.
        double timeStep;

        // Constraint distance vector. Since there are no constraints the
        // vector can be left unassigned.
        std::vector<double> constraintDistances;

        // Current area mismatch.
        double mismatch = computeMismatch(mesh, targetArea);

        // Push current distance from constraint violation into vector.
        constraintDistances.push_back(meshArea*maxMismatch - mismatch);

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

        // Compute the element area fractions.
        boundary.computeAreaFractions();

        // Compute the boundary point normal vectors.
        boundary.computeNormalVectors();

        // Increment the time.
        time += timeStep;

        // Check if the next sample time has been reached.
        while (time >= nextSample)
        {
            // Compute the weighted boundary perimeter.
            double length = computePerimeter(boundary.points, halfHeight, reduce);

            // Record the time, length, and mismatch area.
            times.push_back(time);
            lengths.push_back(length);
            mismatches.push_back(mismatch);

            // Update the time of the next sample.
            nextSample += sampleInterval;

            // Print statistics.
            printf("%6.1f %8.1f %10.4f\n", time, length, mismatch / meshArea);

            // Write level set and boundary segments to file.
            io.saveLevelSetVTK(times.size(), levelSet);
            io.saveBoundarySegmentsTXT(times.size(), boundary);
        }
    }

    // Print results to file.
    FILE *pFile;
    std::ostringstream fileName, num;
    num.str("");
    fileName.str("");
    num << std::fixed << std::setprecision(4) << temperature;
    fileName << "dumbbell_" << num.str() << ".txt";

    pFile = fopen(fileName.str().c_str(), "w");
    for (unsigned int i=0;i<times.size();i++)
        fprintf(pFile, "%lf %lf %lf\n", times[i] - times[0], lengths[i], mismatches[i] / meshArea);
    fclose(pFile);

    std::cout << "\nDone!\n";

    return (EXIT_SUCCESS);
}

// Sensitivity function definition.
double computeSensitivity(const lsm::Coord& coord,
    const lsm::Mesh& mesh, const lsm::LevelSet& levelSet)
{
    /* Interpolate nodal signed distance mismatch to a boundary point
       using inverse squared distance weighting. We are only concerned with
       the sign of the mismatch, i.e. the direction that the boundary should
       move (out or in) in order to reduce the mismatch.
     */

    // Zero the mismatch.
    double mismatch = 0;

    // Find the node that is cloest to the boundary point.
    unsigned int node = mesh.getClosestNode(coord);

    // Loop over node and all of its neighbours.
    for (int i=-1;i<4;i++)
    {
        // Index of neighbour.
        unsigned int n;

        // First test the node itself.
        if (i < 0) n = node;

        // Then its neighbours.
        else n = mesh.nodes[node].neighbours[i];

        // Distance from the boundary point to the node in x & y direction.
        double dx = mesh.nodes[n].coord.x - coord.x;
        double dy = mesh.nodes[n].coord.y - coord.y;

        // Squared distance.
        double rSqd = dx*dx + dy*dy;

        // If boundary point lies exactly on a node then use the sign of
        // the mismatch at that node.
        if (rSqd < 1e-6)
        {
            if (levelSet.target[n] < levelSet.signedDistance[n]) return -1.0;
            else return 1.0;
        }

        // Otherwise update the interpolation estimate.
        else
        {
            if (levelSet.target[n] < levelSet.signedDistance[n]) mismatch -= 1.0 / rSqd;
            else mismatch += 1.0 / rSqd;
        }
    }

    // Return the sign of the interpolated mismatch.
    if (mismatch < 0) return -1.0;
    else return 1.0;
}

// Mismatch function definition.
double computeMismatch(const lsm::Mesh& mesh, const std::vector<double>& targetArea)
{
    double areaMismatch = 0;

    // Compute the total absolute area mismatch.
    for (unsigned int i=0;i<mesh.nElements;i++)
        areaMismatch += std::abs(targetArea[i] - mesh.elements[i].area);

    return areaMismatch;
}

// Perimeter function definition.
double computePerimeter(const std::vector<lsm::BoundaryPoint>& points,
    double halfHeight, double reduce)
{
    double length = 0;

    // Compute the total boundary perimeter.
    for (unsigned int i=0;i<points.size();i++)
    {
        // Reduce the perimeter contribution in the lower half.
        if (points[i].coord.y < halfHeight) length += reduce*points[i].length;
        else length += points[i].length;
    }

    return length;
}
