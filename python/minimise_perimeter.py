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

""" minimise_perimeter.py

    An example code showing a hole shrinking during unconstrained
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
"""
import math
import pyslsm

# Maximum displacement per iteration, in units of the mesh spacing.
# This is the CFL limit.
moveLimit = 0.05

# Set maximum running time.
maxTime = 3000

# Set sampling interval.
sampleInterval = 30

# Set time of the next sample.
nextSample = 30

# Create a hole at position (100, 100) with a radius of 80 grid units.
holes = pyslsm.VectorHole()
holes.append(pyslsm.Hole(100, 100, 80))

# Initialise a 200x200 level set domain.
levelSet = pyslsm.LevelSet(200, 200, holes, moveLimit, 6, True)

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

# Number of cycles since signed distance reinitialisation.
nReinit = 0

# Running time.
runningTime = 0

# Time measurements.
times = pyslsm.VectorDouble()

# Boundary length measurements.
lengths = pyslsm.VectorDouble()

# Boundary curvature measurements.
curvatures = pyslsm.VectorDouble()

# Lambda values for the optimiser.
lambdas = pyslsm.VectorDouble([0])

print("\nStarting unconstrained perimeter minimisation demo...\n")

# Print output header.
print("--------------------------")
print("%6s %8s %10s" % ("Time", "Length", "Curvature"))
print("--------------------------")

# Integrate until we exceed the maximum time.
while runningTime < maxTime:

    # Zero the curvature.
    curvature = 0

    # Initialise the sensitivity object.
    sensitivity = pyslsm.Sensitivity()

    # Initialise the sensitivity callback function.
    cb = pyslsm.Callback()
    cb.callback = boundary.computePerimeter

    # Assign boundary point sensitivities.
    for i in range(0, len(boundary.points)):
        boundary.points[i].sensitivities[0] \
            = sensitivity.computeSensitivity(boundary.points[i], cb.callback)
        curvature += boundary.points[i].sensitivities[0]

    # Compute mean curvature.
    curvature /= len(boundary.points)

    # Time step associated with the iteration.
    timeStep = pyslsm.MutableFloat()

    # Initialise the optimisation object.
    #  Since there are no constraints we pass an empty vector for the
    #  constraint distances argument.
    optimise = pyslsm.Optimise(boundary.points, pyslsm.VectorDouble(), \
        lambdas, timeStep, moveLimit)

    # Perform the optimisation.
    optimise.solve()

    # Extend boundary point velocities to all narrow band nodes.
    levelSet.computeVelocities(boundary.points)

    # Compute gradient of the signed distance function within the narrow band.
    levelSet.computeGradients()

    # Update the level set function.
    isReinitialised = levelSet.update(timeStep.value)

    # Reinitialise the signed distance function, if necessary.
    if (not isReinitialised):
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

    # Compute the boundary point normal vectors.
    boundary.computeNormalVectors(levelSet)

    # Increment the time.
    runningTime += timeStep.value

    # Check if the next sample time has been reached.
    while runningTime >= nextSample:
        # Record the time and boundary length.
        times.append(runningTime)
        lengths.append(boundary.length)
        curvatures.append(curvature)

        # Update the time of the next sample.
        nextSample += sampleInterval

        # Print statistics.
        print("%6.1f %8.1f %10.4f" % (runningTime, boundary.length, curvature))

        # Write level set and boundary segments to file.
        io.saveLevelSetVTK(len(times), levelSet)
        io.saveBoundarySegmentsTXT(len(times), boundary)

# Distance measurements.
distances = pyslsm.VectorDouble()

# Compute the distance moved at each time interval.
for i in range(0, len(times)):
    distances.append((lengths[0] - lengths[i]) / (2 * math.pi))

# Print results to file (distance vs time).
file = open("minimise_perimeter.txt", "w")
for i in range(1, len(times)):
    # Distance and time increments.
    deltaDist = distances[i] - distances[i-1]
    deltaTime = times[i] - times[i-1]
    file.write("%lf %lf %lf %lf %lf\n" % (times[i] - times[0], \
        distances[i], deltaDist / deltaTime, curvatures[i], ((2 * math.pi) / lengths[i])))
file.close()

print("\nDone!")
