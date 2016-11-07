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

#include "slsm.h"

/*! \file shape_match.cpp

    \brief An example showing shape matching optimisation.

    The output file, "shape-match_*.txt", contains the measured mismatch vs
    time data for the optmisation run. Level set information for each sample
    interval is written to ParaView readable VTK files, "level-set_*.vtk".
    Boundary segment data is written to "boundary-segments_*.txt".
 */

// Sensitivity function prototype.
double computeSensitivity(const slsm::Coord&, const slsm::LevelSet&);

// Objective function prototype.
double computeObjective(const slsm::Mesh&, const std::vector<double>&);

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
    double moveLimit = 0.1;

    // Default temperature of the thermal bath.
    double temperature = 0;

    // Override temperature if command-line argument is passed.
    if (argc == 2) temperature = atof(argv[1]);

    // Set maximum running time.
    double maxTime = 100;

    // Set sampling interval.
    double sampleInterval = 1;

    // Set time of the next sample.
    double nextSample = 1;

    // Create a hole at position (100, 100) with a radius of 80 grid units.
    std::vector<slsm::Hole> holes;
    holes.push_back(slsm::Hole(200, 200, 100));

    // Initialise the points vector.
    std::vector<slsm::Coord> points;

    // Read the shape file (assume we're in the root folder).
    std::ifstream shapeFile;
    shapeFile.open("demos/shapes/stanford-bunny.txt");

    if (shapeFile.good())
    {
        // Point coordinates.
        double x, y;

        // Push coordinates into points vector.
        while (shapeFile >> x >> y)
            points.push_back(slsm::Coord({x, y}));
    }
    else
    {
        // Try alternative location (inside demos folder).
        shapeFile.clear();
        shapeFile.open("shapes/stanford-bunny.txt");

        if (shapeFile.good())
        {
            // Point coordinates.
            double x, y;

            // Push coordinates into points vector.
            while (shapeFile >> x >> y)
                points.push_back(slsm::Coord({x, y}));
        }

        else
        {
            std::cerr << "[ERROR]: Invalid shape file!\n";
            exit(EXIT_FAILURE);
        }
    }

    // Close shape file.
    shapeFile.close();

    // Initialise the level set domain.
    slsm::LevelSet levelSet(400, 400, holes, points, moveLimit, 6, true);

    // Initialise io object.
    slsm::InputOutput io;

    // Reinitialise the level set to a signed distance function.
    levelSet.reinitialise();

    // Initialise the boundary object.
    slsm::Boundary boundary;

    // Initialise target area fraction vector.
    std::vector<double> targetArea(levelSet.mesh.nElements);

    // Discretise the target structure.
    boundary.discretise(levelSet, true);

    // Compute the element area fractions.
    levelSet.computeAreaFractions(boundary);

    // Store the target area fractions.
    for (unsigned int i=0;i<levelSet.mesh.nElements;i++)
        targetArea[i] = levelSet.mesh.elements[i].area;

    // Perform initial boundary discretisation.
    boundary.discretise(levelSet);

    // Initialise random number generator.
    slsm::MersenneTwister rng;

    // Number of cycles since signed distance reinitialisation.
    unsigned int nReinit = 0;

    // Running time.
    double time = 0;

    // Time measurements.
    std::vector<double> times;

    // Area mismatch measurements.
    std::vector<double> mismatches;

    /* Lambda values for the optimiser.
       These are reused, i.e. the solution from the current iteration is
       used as an estimate for the next, hence we declare the vector
       outside of the main loop.
     */
    std::vector<double> lambdas(1);

    std::cout << "\nStarting shape matching demo...\n\n";

    // Print output header.
    printf("------------------\n");
    printf("%6s %11s\n", "Time", "Mismatch");
    printf("------------------\n");

    // Integrate until we exceed the maximum time.
    while (time < maxTime)
    {
        // Assign boundary point sensitivities.
        for (unsigned int i=0;i<boundary.points.size();i++)
            boundary.points[i].sensitivities[0] = computeSensitivity(boundary.points[i].coord, levelSet);

        // Initialise the sensitivity object.
        slsm::Sensitivity sensitivity;

        // Apply deterministic Ito correction.
        sensitivity.itoCorrection(boundary, temperature);

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
        boundary.discretise(levelSet);

        // Increment the time.
        time += timeStep;

        // Check if the next sample time has been reached.
        while (time >= nextSample)
        {
            // Compute the element area fractions.
            levelSet.computeAreaFractions(boundary);

            // Current area mismatch.
            double mismatch = computeObjective(levelSet.mesh, targetArea)
                            / (levelSet.mesh.width * levelSet.mesh.height);

            // Record the time and area mismatch.
            times.push_back(time);
            mismatches.push_back(mismatch);

            // Update the time of the next sample.
            nextSample += sampleInterval;

            // Print statistics.
            printf("%6.1f %11.4f\n", time, mismatch);

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
    fileName << "shape-match_" << num.str() << ".txt";

    pFile = fopen(fileName.str().c_str(), "w");
    for (unsigned int i=0;i<times.size();i++)
        fprintf(pFile, "%lf %lf\n", times[i] - times[0], mismatches[i]);
    fclose(pFile);

    std::cout << "\nDone!\n";

    return (EXIT_SUCCESS);
}

// Sensitivity function definition.
double computeSensitivity(const slsm::Coord& coord, const slsm::LevelSet& levelSet)
{
    /* Interpolate nodal signed distance mismatch to a boundary point using
       inverse squared distance weighting. On length scales larger than a
       grid spacing we are only concerned with the sign of the mismatch,
       i.e. the direction that the boundary should move (out or in) in order
       to reduce the mismatch.
     */

    // Zero the mismatch.
    double mismatch = 0;

    // Zero the weighting factor.
    double weight = 0;

    // Find the node that is cloest to the boundary point.
    unsigned int node = levelSet.mesh.getClosestNode(coord);

    // Loop over node and all of its neighbours.
    for (int i=-1;i<4;i++)
    {
        // Index of neighbour.
        unsigned int n;

        // First test the node itself.
        if (i < 0) n = node;

        // Then its neighbours.
        else n = levelSet.mesh.nodes[node].neighbours[i];

        // Distance from the boundary point to the node in x & y direction.
        double dx = levelSet.mesh.nodes[n].coord.x - coord.x;
        double dy = levelSet.mesh.nodes[n].coord.y - coord.y;

        // Squared distance.
        double rSqd = dx*dx + dy*dy;

        // If boundary point lies exactly on a node then use the sign of
        // the mismatch at that node.
        if (rSqd < 1e-6)
        {
            // Calculate nodal mismatch.
            double m = levelSet.target[n] - levelSet.signedDistance[n];

            // Smooth mismatch over a grid spacing.
            if (std::abs(m) < 1.0) return m;

            // Return the sign of the mismatch.
            if (m < 0) return -1.0;
            else return 1.0;
        }

        // Otherwise, update the interpolation estimate.
        else
        {
            mismatch += (levelSet.target[n] - levelSet.signedDistance[n]) / rSqd;
            weight   += 1.0 / rSqd;
        }
    }

    // Compute weighted mismatch.
    mismatch /= weight;

    // Smooth mismatch over a grid spacing.
    if (std::abs(mismatch) < 1.0) return mismatch;

    // Return the sign of the interpolated mismatch.
    if (mismatch < 0) return -1.0;
    else return 1.0;
}

// Objective function definition.
double computeObjective(const slsm::Mesh& mesh, const std::vector<double>& targetArea)
{
    double areaMismatch = 0;

    // Compute the total absolute area mismatch.
    for (unsigned int i=0;i<mesh.nElements;i++)
        areaMismatch += std::abs(targetArea[i] - mesh.elements[i].area);

    return areaMismatch;
}
