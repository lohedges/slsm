/*
  Copyright (c) 2015-2017 Lester Hedges <lester.hedges+lsm@gmail.com>

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

#include "slsm.h"

#ifndef M_PI
    #define M_PI 3.1415926535897932384626433832795
#endif

/*! \file bimodal_brolly.cpp

    \brief An umbrella sampling demo for the bimodal objective function.

    Here we construct a simple example system with two minima separated by
    a small free energy barrier. The matched shape is a dumbbell constructed
    from two horizontally offset, overlapping circles. The position of the
    right-hand dumbbell lobe is shifted vertically by a small amount relative
    to the left.

    We construct two minima by making the perimeter (objective) sensitivities
    a function of the y coordinate of each boundary point, mimicking a
    gravitational field. Sensitivities are linearly reduced between top and
    bottom of the left-hand lobe, hence within the left lobe it's possible to
    form a circle with a smaller perimeter at the same cost.

    Umbrella sampling allows us to sample equilibrium states that are low
    in probability (high in free energy) by constraining the system using
    a harmonic bias potential. By combining sampling data from different
    umbrella windows we can calculate the free energy profile for the
    transition between the left and right lobe. Our bias potential constrains
    the x centre of mass of the shape <x>, i.e. the bias is (k/2)*(x_s - x_i)^2,
    where k is the spring constant, x_s is the <x> for the current sample,
    and x_i is the <x> for umbrella window i.

    The output file, "brolly_*.txt", contains the measured x centre of mass,
    perimeter and mismatch vs time data for the umbrella sampling run. Level set
    information for each sample interval is written as a binary restart file,
    "level-set_*.bin".  Boundary segment data is written to
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
double computePerimeterWeight(double y);

// Calculate the "centre of mass" of the boundary".
void computeCentreOfMass(const std::vector<slsm::BoundaryPoint>&, double&, double&);

// Bias potential function prototype.
double computeBiasPotential(double, double, double);

// Acceptance function prototype.
bool isAccepted(double, double, double, slsm::MersenneTwister&);

// Parse arguments from the command-line.
void parseCommandLineArguments(int, char**, double&, double&,
    double&, double&, double&, unsigned int&, unsigned int&, char*, bool&);

// Print help message to stdout.
void printHelpMessage();

// GLOBALS
unsigned int nDiscrete = 10;                // Boundary integral discretisation factor.
double uppergravityCutOff;                  // The maximum y coordinate at which gravity is active.
double lowergravityCutOff;                  // The minimum y coordinate at which gravity is active.
double gravityRange;                        // The vertical separation between dumbbell lobes.
double reduce;                              // Sensitivity reduction factor.
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

    // Set default paramaters.

    // Temperature of the thermal bath.
    double temperature = 0.05;

    // Gravitational multiplier (must be between zero and 1).
    double gravityMult = 0.02;

    // Centre of the bias potential.
    double centre = 10;

    // Spring constant for harmonic bias.
    double spring = 0.2;

    // Time interval for umbrella sampling.
    double umbrellaInterval = 1;

    // Number of umbrella sampling steps per sample.
    unsigned int sampleInterval = 10;

    // Number of samples.
    unsigned int nSamples = 10000;

    // File name for the starting configuration.
    char restart[100];

    // Whether a restart file is being used.
    bool isRestart = false;

	// Read command-line arguments.
    parseCommandLineArguments(argc, argv, temperature, gravityMult, centre,
        spring, umbrellaInterval, sampleInterval, nSamples, restart, isRestart);

    // Print parameters.
    printf("\nParameters:\n");
    printf(" temperature       = %3.2f\n", temperature);
    printf(" gravityMult       = %3.2f\n", gravityMult);
    printf(" centre            = %3.2f\n", centre);
    printf(" spring            = %3.2f\n", spring);
    printf(" umbrellaInterval  = %3.2f\n", umbrellaInterval);
    printf(" sampleInterval    = %d\n", sampleInterval);
    printf(" nSamples          = %d\n", nSamples);
    if (isRestart)
        printf(" restart       = %s\n", restart);

    // Inverse temperature.
    double beta = 1.0 / temperature;

    // Width and height of the mesh.
    double width  = 70;
    double height = 60;

    // Hole vectors.
    std::vector<slsm::Hole> initialHoles;
    std::vector<slsm::Hole> targetHoles;

    // Create a dumbbell from two horizontally offset holes.

    // Centre the dumbell horizontally and vertically.
    double xCentre = 0.5*width;
    double yCentre = 0.5*height;

    // Vertical offset between dumbbell lobes.
    double yOffset = 1;

    // Horizontal separation of lobe centres is 2.xOffset.
    double xOffset = 10;

    // Lobe radius.
    double radius  = 20;

    // Push holes into vectors.
    targetHoles.push_back(slsm::Hole(xCentre+xOffset, yCentre+yOffset, radius));
    targetHoles.push_back(slsm::Hole(xCentre-xOffset, yCentre, radius));

    double aa = xOffset / radius;
    double area = (M_PI/2.0) + asin(aa) + aa*sqrt(1-aa*aa);
    area *= 2.0 * radius * radius;

    // Set maximumum area mismatch (as fraction of total area).
    double minSizeOfShape = 0.2 / 3.0;
    minSizeOfShape *= width*height/area;

    // Set the mismatch constraint.
    double maxMismatch = area/width/height * (1.0 - minSizeOfShape);

    // Store dumbbell data (needed for sensitivity callback).
    // They are used for implementing a cutoff for the "gravity".
    uppergravityCutOff = targetHoles[0].coord.y + radius;
    lowergravityCutOff = targetHoles[0].coord.y - radius;
    gravityRange = uppergravityCutOff - lowergravityCutOff;

    /* Physics stuff:
        The "energy" of a boundary segment is l.g.y where l is length, g is gravity, y is height.
        The original idea is that g.y=1      at uppergravityCutOff
                             and  g.y=reduce at lowergravityCutOff
        with linear interpolation in between
          hence g.(upper-lower) = 1-reduce
                g = (1-reduce)/(upper-lower) is the g-field
          or
                reduce = 1 - g.(upper-lower), which must be positive
          hence
                0 < g < 1/(upper-lower)
     */

    // The variable gravityMult is g.(upper-lower).
    reduce = 1.0 - gravityMult;

    // Position and size of the initial trial circle.
    double initialHoleX = xCentre;
    double initialHoleY = 0.5*yCentre;
    double initialHoleRad = 0.5*radius;

    // Push the trial shape into the vector.
    initialHoles.push_back(slsm::Hole(initialHoleX, initialHoleY, initialHoleRad));

    // Initialise the level set domain.
    slsm::LevelSet levelSet(width, height, initialHoles, targetHoles, moveLimit, 6, true);

    // Store the mesh area.
    double meshArea = levelSet.mesh.width * levelSet.mesh.height;

    // Initialise io object.
    slsm::InputOutput io;

    // Read binary restart file.
    if (isRestart)
    {
        std::ostringstream restartFile(restart);
        io.loadLevelSetBIN(restartFile, levelSet);
    }

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
    double runningTime = 0;

    // Backup the signed distance function.
    std::vector<double> signedDistance = levelSet.signedDistance;

    // Compute the initial centre of mass.
    double xCentreOfMass, tmp;
    computeCentreOfMass(boundary.points, xCentreOfMass, tmp);

    // Compute the initial bias potential.
    double biasPotential = computeBiasPotential(xCentreOfMass, centre, spring);

    /* Lambda values for the optimiser.
       These are reused, i.e. the solution from the current iteration is
       used as an estimate for the next, hence we declare the vector
       outside of the main loop.
     */
    std::vector<double> lambdas(2);

    // Set up file name.
    std::ostringstream fileName, fileName2, fileName3;
    fileName.precision(2);
    fileName2.precision(2);
    fileName3.precision(2);
    fileName  << std::fixed;
    fileName2 << std::fixed;
    fileName3 << std::fixed;
    fileName  << "brolly_" << centre << ".txt";
    fileName2 << "level-set_" << centre << ".bin";
    fileName3 << "boundary-segments_" << centre << ".txt";

    // Wipe existing log file.
    FILE *pFile;
    pFile = fopen(fileName.str().c_str(), "w");
    fclose(pFile);

    // Number of accepted trials and total trials.
    unsigned int nAccept = 0;
    unsigned int nTrials = 0;

    std::cout << "\nStarting umbrella sampling demo...\n\n";

    // Print output header.
    printf("----------------------------------------------------\n");
    printf("%8s %10s %10s %10s %10s\n", "Time", "<x>", "Length", "Mismatch", "Accept");
    printf("----------------------------------------------------\n");

    for (unsigned int i=0;i<nSamples;i++)
    {
        for (unsigned int j=0;j<sampleInterval;j++)
        {
            // Zero the sample interval time.
            double sampleTime = 0;

            // Integrate until we exceed umbrella sampling interval.
            while (sampleTime < umbrellaInterval)
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

                // Set constraint type.
                std::vector<bool> isEquality;
                isEquality.push_back(false);

                /* Initialise the optimisation object.

                The Optimise class is a lightweight object so there is no cost for
                reinitialising at every iteration. A smart compiler will optimise
                this anyway, i.e. the same memory space will be reused. It is better
                to place objects in the correct scope in order to aid readability
                and to avoid unintended name clashes, etc.
                */
                slsm::Optimise optimise(boundary.points, constraintDistances,
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

                // Compute the element area fractions.
                levelSet.computeAreaFractions(boundary);

                // Compute the boundary point normal vectors.
                boundary.computeNormalVectors(levelSet);

                // Increment the time.
                sampleTime += timeStep;
            }

            double xPos, yPos;

            // Compute the centre of mass.
            computeCentreOfMass(boundary.points, xPos, yPos);
            xPos -= 0.5*width;

            // Compute the trial bias potential.
            double biasPotentialTrial = computeBiasPotential(xPos, centre, spring);

            // See if sample is accepted.
            if (isAccepted(biasPotentialTrial, biasPotential, beta, rng))
            {
                // Store updated measurements.
                xCentreOfMass = xPos;
                biasPotential = biasPotentialTrial;

                // Backup the current signed distance function.
                signedDistance = levelSet.signedDistance;

                nAccept++;
            }
            else
            {
                // Reset the signed distance function.
                levelSet.signedDistance = signedDistance;

                // Reinitialise the signed distance function.
                levelSet.reinitialise();
                nReinit = 0;

                // Compute the new discretised boundary.
                boundary.discretise(levelSet);

                // Compute the element area fractions.
                levelSet.computeAreaFractions(boundary);

                // Compute the boundary point normal vectors.
                boundary.computeNormalVectors(levelSet);
            }

            // Update the total running time.
            runningTime += sampleTime;

            // Increment the number of trials.
            nTrials++;
        }

        // Current area mismatch.
        double mismatch = computeMismatch(levelSet.mesh, targetArea);

        // Current weighted perimeter.
        double length = computePerimeter(boundary.points);

        // Print sample to stdout.
        printf("%6.2e %10.4f %10.4f %10.4f %10.4f\n",
            runningTime, xCentreOfMass, length, mismatch / meshArea, ((double) nAccept / nTrials));

        // Write sample to file.
        pFile = fopen(fileName.str().c_str(), "a");
        fprintf(pFile, "%e %lf %lf %lf %lf\n",
            runningTime, xCentreOfMass, length, mismatch / meshArea, ((double) nAccept / nTrials));
        fclose(pFile);

        // Write level set and boundary segments to file.
        io.saveLevelSetBIN(fileName2, levelSet);
        io.saveBoundarySegmentsTXT(fileName3, boundary);
    }

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
    if      (y > uppergravityCutOff) return 1.0;
    else if (y < lowergravityCutOff) return reduce;
    else
    {
        // Fractional distance along the gravity range.
        double dy = (uppergravityCutOff - y) / gravityRange;

        return (1.0 - dy*(1.0 - reduce));
    }
}

// Boundary "centre of mass" function definition.
void computeCentreOfMass(const std::vector<slsm::BoundaryPoint>& points, double& x, double& y)
{
    // Zero variables.
    double length = 0;
    x = 0;
    y = 0;

    // Compute the "centre of mass" of the boundary.
    for (unsigned int i=0;i<points.size();i++)
    {
        // Compute the weighted point length.
        double ll = computePointLength(points[i]);

        // Update variables.
        length += 0.5*ll;
        x += points[i].coord.x * 0.5*ll;
        y += points[i].coord.y * 0.5*ll;
    }

    // Normalise.
    x /= length;
    y /= length;
}

// Bias potential function definition.
double computeBiasPotential(double value, double centre, double spring)
{
    return 0.5*spring*(value - centre)*(value - centre);
}

// Acceptance function definition.
bool isAccepted(double currentbias, double previousBias, double beta, slsm::MersenneTwister& rng)
{
    if (rng() < exp(-beta*(currentbias - previousBias))) return true;
    else return false;
}

// Parse arguments from the command-line.
void parseCommandLineArguments(int argc, char **argv, double& temperature,
    double& gravityMult, double& centre, double& spring, double& umbrellaInterval,
    unsigned int& sampleInterval, unsigned int& nSamples, char* restart, bool& isRestart)
{
    unsigned int i = 1;

    if (argc == 2)
    {
        if (strcmp(argv[1],"-h") == 0 || strcmp(argv[1],"--help") == 0)
        {
            printHelpMessage();
            exit(EXIT_SUCCESS);
        }
    }

    else
    {
        while (i < argc)
        {
            if (strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--temperature") == 0)
            {
                i++;
                temperature = atof(argv[i]);
            }

            else if (strcmp(argv[i], "-g") == 0 || strcmp(argv[i], "--gravity") == 0)
            {
                i++;
                gravityMult = atof(argv[i]);
            }

            else if (strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "--centre") == 0)
            {
                i++;
                centre = atof(argv[i]);
            }

            else if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--spring") == 0)
            {
                i++;
                spring = atof(argv[i]);
            }

            else if (strcmp(argv[i], "-ui") == 0 || strcmp(argv[i], "--unmbrella-interval") == 0)
            {
                i++;
                umbrellaInterval = atof(argv[i]);
            }

            else if (strcmp(argv[i], "-si") == 0 || strcmp(argv[i], "--sample-interval") == 0)
            {
                i++;
                sampleInterval = atoi(argv[i]);
            }

            else if (strcmp(argv[i], "-n") == 0 || strcmp(argv[i], "--number-samples") == 0)
            {
                i++;
                nSamples = atoi(argv[i]);
            }

            else if (strcmp(argv[i], "-r") == 0 || strcmp(argv[i], "--restart") == 0)
            {
                i++;
                isRestart = true;
                strcpy(restart, argv[i]);
            }

            else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0)
            {
                printHelpMessage();
                exit(EXIT_SUCCESS);
            }

            else
            {
                fprintf(stderr, "\n[ERROR]: Unknown command-line option: %s\n", argv[i]);
                fprintf(stderr, "For help run \"bimodal_brolly -h\"\n");
                exit(EXIT_FAILURE);
            }

            i++;
        }
    }
}

// Print help message to stdout.
void printHelpMessage()
{
    puts("\nbimodal_brolly\n\n"

         "Available options:\n"
         " -h/--help                          : Print this help information\n"
         " -t/--temperature <double>          : Temperature of the thermal bath\n"
         " -g/--gravity <double>              : Gravitational multiplier\n"
         " -c/--centre <double>               : Centre of the umbrella potential\n"
         " -s/--spring <double>               : Spring constant\n"
         " -ui/--umbrella-interval <double>   : Length of trial trajectory\n"
         " -si/--sample-interval <int>        : Sampling frequency\n"
         " -n/--number-samples <int>          : Total number of samples\n"
         " -r/--restart <string>              : Location of restart file\n");
}
