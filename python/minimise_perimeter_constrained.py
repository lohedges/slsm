#  Copyright (c) 2015-2017 Lester Hedges <lester.hedges+slsm@gmail.com>
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program. If not, see <http://www.gnu.org/licenses/>.

""" minimise_perimeter_constrained.py

    An example code showing perimeter minimisation with noise.
    Minimisation is performed subject to a lower bound on the area (hole area).
    (Since the shape is cut from a filled domain this corresponds to an upper
    bound on the area outside of the shape.)

    The output file, "perimeter_*.txt", contains the measured boundary length vs
    time data for the optmisation run. Level set information for each sample
    interval is written to ParaView readable VTK files, "level-set_*.vtk".
    Boundary segment data is written to "boundary-segments_*.txt".
"""

import math
import pyslsm
import sys

# Maximum displacement per iteration, in units of the mesh spacing.
# This is the CFL limit.
moveLimit = 0.1

# Default temperature of the thermal bath.
temperature = 0.1

# Override temperature if command-line argument is passed.
if (len(sys.argv) == 2):
    temperature = float(sys.argv[1])

# Set maximum running time.
maxTime = 50

# Set maximum material area.
maxArea = 0.6

# Set sampling interval.
sampleInterval = 0.5

# Set time of the next sample.
nextSample = 0.5

# Initialise a 200x200 level set domain.
levelSet = pyslsm.LevelSet(200, 200, moveLimit, 6, True)

# Calculate the area of the mesh.
meshArea = levelSet.mesh.width * levelSet.mesh.height

# Create a solid slab of material with a small square in the middle.
for i in range(0, levelSet.mesh.nNodes):
    x = math.floor(levelSet.mesh.nodes[i].coord.x)
    y = math.floor(levelSet.mesh.nodes[i].coord.y)

    # Cut out square hole.
    if (x >= 90 and x <= 110 and y >=90 and y <= 110):
        levelSet.signedDistance[i] = -1
    else:
        levelSet.signedDistance[i] = 1

# Initialise io object.
io = pyslsm.InputOutput()

# Reinitialise the level set to a signed distance function.
levelSet.reinitialise()

# Initialise the boundary object.
boundary = pyslsm.Boundary()

# Perform initial boundary discretisation.
boundary.discretise(levelSet)

# Compute the initial boundary point normal vectors.
boundary.computeNormalVectors(levelSet)

# Initialise random number generator.
rng = pyslsm.MersenneTwister()

# Number of cycles since signed distance reinitialisation.
nReinit = 0

# Running time.
runningTime = 0

# Time measurements.
times = pyslsm.VectorDouble()

# Boundary length measurements (objective).
lengths = pyslsm.VectorDouble()

# Material area fraction measurements (constraint).
areas = pyslsm.VectorDouble()

# Lambda values for the optimiser.
#  These are reused, i.e. the solution from the current iteration is
#  used as an estimate for the next, hence we declare the vector
#  outside of the main loop.
lambdas = pyslsm.VectorDouble([0, 0])

print("\nStarting constrained perimeter minimisation demo...\n")

# Print output header.
print("-------------------------")
print("%9s %8s %6s" % ("Time", "Length", "Area"))
print("-------------------------")

# Integrate until we exceed the maximum time.
while runningTime < maxTime:

    # Initialise the sensitivity object.
    sensitivity = pyslsm.Sensitivity()

    # Initialise the sensitivity callback function.
    cb = pyslsm.Callback()
    cb.callback = boundary.computePerimeter

    # Assign boundary point sensitivities.
    for i in range(0, len(boundary.points)):
        boundary.points[i].sensitivities[0] \
            = sensitivity.computeSensitivity(boundary.points[i], cb.callback)
        boundary.points[i].sensitivities[1] = -1.0

    # Apply deterministic Ito correction.
    sensitivity.itoCorrection(boundary, temperature)

    # Time step associated with the iteration.
    timeStep = pyslsm.MutableFloat()

    # Constraint distance vector.
    constraintDistances = pyslsm.VectorDouble()

    # Push current distance from constraint violation into vector.
    constraintDistances.append(meshArea*maxArea - levelSet.area)

    # Initialise the optimisation object.
    #  Since there are no constraints we pass an empty vector for the
    #  constraint distances argument.
    optimise = pyslsm.Optimise(boundary.points, constraintDistances, \
        lambdas, timeStep, moveLimit)

    # Perform the optimisation.
    optimise.solve()

    # Extend boundary point velocities to all narrow band nodes.
    levelSet.computeVelocities(boundary.points, timeStep, temperature, rng)

    # Compute gradient of the signed distance function within the narrow band.
    levelSet.computeGradients()

    # Update the level set function.
    isReinitialised = levelSet.update(timeStep.value)

    # Reinitialise the signed distance function, if necessary.
    if not isReinitialised:
        # Reinitialise at least every 20 iterations.
        if (nReinit == 20):
            levelSet.reinitialise()
            nReinit = 0
    else:
        nReinit = 0

    # Increment the number of steps since reinitialisation.
    nReinit += 1

    # Compute the new discretised boundary.
    boundary.discretise(levelSet)

    # Compute the element area fractions.
    levelSet.computeAreaFractions(boundary)

    # Compute the boundary point normal vectors.
    boundary.computeNormalVectors(levelSet)

    # Increment the time.
    runningTime += timeStep.value

    # Check if the next sample time has been reached.
    while runningTime >= nextSample:
        # Record the time, boundary length, and material area.
        times.append(runningTime)
        lengths.append(boundary.length)
        areas.append(levelSet.area / meshArea)

        # Update the time of the next sample.
        nextSample += sampleInterval

        # Print statistics.
        print("%9.2f %8.1f %6.2f" % (runningTime, boundary.length, levelSet.area / meshArea))

        # Write level set and boundary segments to file.
        io.saveLevelSetVTK(len(times), levelSet)
        io.saveBoundarySegmentsTXT(len(times), boundary)

# Print results to file (distance vs time).
fileName = "perimeter_%5.4f.txt" % (temperature)
file = open(fileName, "w")
for i in range(0, len(times)):
    file.write("%lf %lf %lf\n" % (times[i] - times[0], lengths[i], areas[i]))
file.close()

print("\nDone!")
