/*
  Copyright (c) 2015-2017 Lester Hedges <lester.hedges+slsm@gmail.com>

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

/*! \file dumbbell.cpp

    \brief An example showing perimeter minimisation with a shape matching
    constraint.

    Here we construct a simple example system with two minima separated by
    a free energy barrier. The matched shape is a narrow-necked dumbbell
    constructed from two vertically offset, overlapping circles. The initial
    configuration is a circle centred in the upper lobe of the dumbbell.
    We construct two minima by making the perimeter (objective) sensitivities
    a function of the y coordinate of each boundary point. Sensitivities
    are linearly reduced between the centres of the upper and lower lobes,
    hence within the lowesr lobe it's possible to form a circle with a smaller
    perimeter at the same cost.

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

// FUNCTION PROTOTYPES

// Constraint sensitivity function prototype.
double computeConstraintSensitivity(const slsm::Coord&, const slsm::LevelSet&);

// Mismatch function prototype.
double computeMismatch(const slsm::Mesh&, const std::vector<double>&);

// Perimeter function prototype.
double computePerimeter(const std::vector<slsm::BoundaryPoint>&);

// Boundary point length function prototype.
double computePointLength(const slsm::BoundaryPoint& point);

// Perimeter weight function prototype.
double computePerimeterWeight(double);

// GLOBALS
unsigned int nDiscrete = 10;                // Boundary integral discretisation factor.
double upperLobeCentre;                     // The y coordinate of the upper dumbbell lobe.
double lowerLobeCentre;                     // The y coordinate of the lower dumbbell lobe.
double lobeSeparation;                      // The vertical separation between dumbbell lobes.
double reduce = 0.65;                       // Sensitivity reduction factor.
std::vector<slsm::BoundaryPoint>* points;   // Pointer to the boundary points vector.

// MAIN FUNCTION

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
    double temperature = 0.02;

    // Override temperature if command-line argument is passed.
    if (argc > 1) temperature = atof(argv[1]);

    // Override sensitivity reduction factor if command-line argument is passed.
    if (argc > 2) reduce = atof(argv[2]);

    // Set maximum running time.
    double maxTime = 8000;

    // Set maximumum area mismatch.
    double maxMismatch = 0.2;

    // Set sampling interval.
    double sampleInterval = 40;

    // Set time of the next sample.
    double nextSample = 40;

    // Hole vectors.
    std::vector<slsm::Hole> initialHoles;
    std::vector<slsm::Hole> targetHoles;

    // Create a dumbbell from two vertically offset holes.
    targetHoles.push_back(slsm::Hole(50, 69, 20));
    targetHoles.push_back(slsm::Hole(50, 31, 20));

    // Store dumbbell data (needed for sensitivity callback).
    upperLobeCentre = targetHoles[0].coord.y;
    lowerLobeCentre = targetHoles[1].coord.y;
    lobeSeparation = upperLobeCentre - lowerLobeCentre;

    // Initialise the system with a circle in the upper lobe.
    initialHoles.push_back(slsm::Hole(50, 69, 15));

    // Initialise a 100x100 level set domain.
    slsm::LevelSet levelSet(100, 100, initialHoles, targetHoles, moveLimit, 6, true);

    // Store the mesh area.
    double meshArea = levelSet.mesh.width * levelSet.mesh.height;

    // Read starting configuration from file.
    if (argc > 3)
    {
        // Open the starting configuration.
        std::ifstream inputFile(argv[3]);

        // Read signed distance from file.
        if (inputFile.good())
        {
            unsigned int i=0;
            while (inputFile >> levelSet.signedDistance[i])
                i++;
        }
        else
        {
            std::cout << "Invalid configuration file!\n";
            exit(EXIT_FAILURE);
        }
    }

    // Initialise io object.
    slsm::InputOutput io;

    // Reinitialise the level set to a signed distance function.
    levelSet.reinitialise();

    // Initialise the boundary object.
    slsm::Boundary boundary;

    // Initialise boundary points pointer.
    points = &boundary.points;

    // Initialise target area fraction vector.
    std::vector<double> targetArea(levelSet.mesh.nElements);

    // Discretise the target structure.
    boundary.discretise(levelSet, true);
    io.saveBoundarySegmentsTXT(0, boundary);

    // Compute the element area fractions.
    levelSet.computeAreaFractions(boundary);

    // Store the target area fractions.
    for (unsigned int i=0;i<levelSet.mesh.nElements;i++)
        targetArea[i] = levelSet.mesh.elements[i].area;

    // Perform initial boundary discretisation.
    boundary.discretise(levelSet);

    // Compute the element area fractions.
    levelSet.computeAreaFractions(boundary);

    // Compute the initial boundary point normal vectors.
    boundary.computeNormalVectors(levelSet);

    // Initialise random number generator.
    slsm::MersenneTwister rng;

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
        slsm::Sensitivity sensitivity;

        // Initialise the sensitivity callback.
        using namespace std::placeholders;
        slsm::SensitivityCallback callback = std::bind(&computePointLength, _1);

        // Assign boundary point sensitivities.
        for (unsigned int i=0;i<boundary.points.size();i++)
        {
            boundary.points[i].sensitivities[0] =
                sensitivity.computeSensitivity(boundary.points[i], callback);
            boundary.points[i].sensitivities[1] =
                computeConstraintSensitivity(boundary.points[i].coord, levelSet);
        }

        // Apply deterministic Ito correction.
        sensitivity.itoCorrection(boundary, temperature);

        // Time step associated with the iteration.
        double timeStep;

        // Constraint distance vector.
        std::vector<double> constraintDistances;

        // Current area mismatch.
        double mismatch = computeMismatch(levelSet.mesh, targetArea);

        // Push current distance from constraint violation into vector.
        constraintDistances.push_back(meshArea*maxMismatch - mismatch);

        // Initialise the optimisation object.
        slsm::Optimise optimise(boundary.points, constraintDistances,
            lambdas, timeStep, levelSet.moveLimit);

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

        // Compute the element area fractions.
        levelSet.computeAreaFractions(boundary);

        // Compute the boundary point normal vectors.
        boundary.computeNormalVectors(levelSet);

        // Increment the time.
        time += timeStep;

        // Check if the next sample time has been reached.
        while (time >= nextSample)
        {
            // Compute the weighted boundary perimeter.
            double length = computePerimeter(boundary.points);

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

// FUNCTION DEFINITIONS

// Constraint sensitivity function definition.
double computeConstraintSensitivity(const slsm::Coord& coord, const slsm::LevelSet& levelSet)
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

// Mismatch function definition.
double computeMismatch(const slsm::Mesh& mesh, const std::vector<double>& targetArea)
{
    double areaMismatch = 0;

    // Compute the total absolute area mismatch.
    for (unsigned int i=0;i<mesh.nElements;i++)
        areaMismatch += std::abs(targetArea[i] - mesh.elements[i].area);

    return areaMismatch;
}

// Perimeter function definition.
double computePerimeter(const std::vector<slsm::BoundaryPoint>& points)
{
    double length = 0;

    // Compute the total weighted boundary length.
    for (unsigned int i=0;i<points.size();i++)
        length += 0.5*computePointLength(points[i]);

    return length;
}

// Boundary point length function definition.
double computePointLength(const slsm::BoundaryPoint& point)
{
    double length = 0;

    // Store coordinates of the point.
    double x1 = point.coord.x;
    double y1 = point.coord.y;

    // Sum the distance to each neighbour.
    for (unsigned int i=0;i<point.nNeighbours;i++)
    {
        // Store coordinates of the neighbouring point.
        double x2 = (*points)[point.neighbours[i]].coord.x;
        double y2 = (*points)[point.neighbours[i]].coord.y;

        // Compute separation components.
        double dx = x2 - x1;
        double dy = y2 - y1;

        // Compute segment length.
        double len = sqrt(dx*dx + dy*dy) / nDiscrete;

        // Perform discrete boundary (line) integral.
        for (unsigned int j=0;j<nDiscrete;j++)
        {
            // Compute y position along segment.
            double y = y1 + (j+0.5)*dy/nDiscrete;

            // Add weighted segment length.
            length += computePerimeterWeight(y)*len;
        }
    }

    return length;
}

// Perimeter weight function definition.
double computePerimeterWeight(double y)
{
    if      (y > upperLobeCentre) return 1.0;
    else if (y < lowerLobeCentre) return reduce;
    else
    {
        // Fractional distance from the centre of the upper dumbbell lobe.
        double dy = (upperLobeCentre - y) / lobeSeparation;

        return (1.0 - dy*(1.0 - reduce));
    }
}
