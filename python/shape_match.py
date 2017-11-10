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

""" shape_match.py

    An example showing shape matching optimisation.

    The output file, "shape-match_*.txt", contains the measured mismatch vs
    time data for the optmisation run. Level set information for each sample
    interval is written to ParaView readable VTK files, "level-set_*.vtk".
    Boundary segment data is written to "boundary-segments_*.txt".
"""

import math
import pyslsm
import sys

# Sensitivity function.
def computeSensitivity(coord, levelSet):
    """
    Interpolate nodal signed distance mismatch to a boundary point using
    inverse squared distance weighting. On length scales larger than a
    grid spacing we are only concerned with the sign of the mismatch,
    i.e. the direction that the boundary should move (out or in) in order
    to reduce the mismatch.
    """

    # Zero the mismatch.
    mismatch = 0

    # Zero the weighting factor.
    weight = 0

    # Find the node that is cloest to the boundary point.
    node = levelSet.mesh.getClosestNode(coord)

    # Loop over node and all of its neighbours.
    for i in range(-1,3):
        # Work out index of neighbour.

        # First test the node itself.
        if i < 0:
            n = node
        # Then its neighbours.
        else:
            n = levelSet.mesh.nodes[node].neighbours[i]

        # Distance from the boundary point to the node in x & y direction.
        dx = levelSet.mesh.nodes[n].coord.x - coord.x
        dy = levelSet.mesh.nodes[n].coord.y - coord.y

        # Squared distance.
        rSqd = dx*dx + dy*dy

        # If boundary point lies exactly on a node then use the sign of
        # the mismatch at that node.
        if rSqd < 1e-6:
            # Calculate nodal mismatch.
            m = levelSet.target[n] - levelSet.signedDistance[n]

            # Smooth mismatch over a grid spacing.
            if abs(m) < 1.0:
                return m

            # Return the sign of the mismatch.
            if m < 0:
             return -1.0

            else:
             return 1.0

        # Otherwise, update the interpolation estimate.
        else:
            mismatch += (levelSet.target[n] - levelSet.signedDistance[n]) / rSqd
            weight   += 1.0 / rSqd

    # Compute weighted mismatch.
    mismatch /= weight

    # Smooth mismatch over a grid spacing.
    if abs(mismatch) < 1.0:
        return mismatch

    # Return the sign of the interpolated mismatch.
    if mismatch < 0:
        return -1.0

    else:
        return 1.0

# Objective function.
def computeObjective(mesh, targetArea):
    """
    Compute the total area mismatch between the level set domain
    and its target value.
    """

    # Zero the mismatch.
    areaMismatch = 0

    # Compute the total absolute area mismatch.
    for i in range(0, mesh.nElements):
        areaMismatch += abs(targetArea[i] - mesh.elements[i].area)

    return areaMismatch

# Maximum displacement per iteration, in units of the mesh spacing.
# This is the CFL limit.
moveLimit = 0.1

# Default temperature of the thermal bath.
temperature = 0

# Override temperature if command-line argument is passed.
if (len(sys.argv) == 2):
    temperature = float(sys.argv[1])

# Set maximum running time.
maxTime = 100

# Set sampling interval.
sampleInterval = 1

# Set time of the next sample.
nextSample = 1

# Create a hole at position (200, 200) with a radius of 100 grid units.
holes = pyslsm.VectorHole()
holes.append(pyslsm.Hole(200, 200, 100))

# Initialise the points vector.
points = pyslsm.VectorCoord()

# Read the shape file (assume we're in the root folder).
try:
    file = open("demos/shapes/stanford-bunny.txt", "r")

    for line in file.readlines():
        cols = line.split()
        points.append(pyslsm.Coord(float(cols[0]), float(cols[1])))

except (OSError, IOError):
    print("Shape file not found!")
    sys.exit()

# Initialise the level set domain.
levelSet = pyslsm.LevelSet(400, 400, holes, points, moveLimit, 6, True)

# Initialise io object.
io = pyslsm.InputOutput()

# Reinitialise the level set to a signed distance function.
levelSet.reinitialise()

# Initialise the boundary object.
boundary = pyslsm.Boundary()

# Initialise target area fraction vector.
targetArea = pyslsm.VectorDouble()

# Discretise the target structure.
boundary.discretise(levelSet, True)

# Compute the element area fractions.
levelSet.computeAreaFractions(boundary)

# Store the target area fractions.
for i in range(0, levelSet.mesh.nElements):
    targetArea.append(levelSet.mesh.elements[i].area)

# Perform initial boundary discretisation.
boundary.discretise(levelSet)

# Initialise random number generator.
rng = pyslsm.MersenneTwister()

# Number of cycles since signed distance reinitialisation.
nReinit = 0

# Running time.
runningTime = 0

# Time measurements.
times = pyslsm.VectorDouble()

# Area mismatch measurements.
mismatches = pyslsm.VectorDouble()

# Lambda values for the optimiser.
lambdas = pyslsm.VectorDouble([0])

print("Starting shape matching demo...\n")

# Print output header.
print("------------------")
print("%6s %11s" % ("Time", "Mismatch"))
print("------------------")

# Integrate until we exceed the maximum time.
while runningTime < maxTime:

    # Assign boundary point sensitivities.
    for i in range(0, len(boundary.points)):
        boundary.points[i].sensitivities[0] = computeSensitivity(boundary.points[i].coord, levelSet)

    # Initialise the sensitivity object.
    sensitivity = pyslsm.Sensitivity()

    # Apply deterministic Ito correction.
    sensitivity.itoCorrection(boundary, temperature)

    # Time step associated with the iteration.
    timeStep = pyslsm.MutableFloat()

    # Initialise the optimisation object.
    optimise = pyslsm.Optimise(boundary.points, pyslsm.VectorDouble(), \
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
        if nReinit == 20:
            levelSet.reinitialise()
            nReinit = 0
    else:
        nReinit = 0

    # Increment the number of steps since reinitialisation.
    nReinit += 1

    # Compute the new discretised boundary.
    boundary.discretise(levelSet)

    # Increment the time.
    runningTime += timeStep.value

    # Check if the next sample time has been reached.
    while runningTime >= nextSample:
        # Compute the element area fractions.
        levelSet.computeAreaFractions(boundary)

        # Current area mismatch.
        mismatch = computeObjective(levelSet.mesh, targetArea) \
                  / (levelSet.mesh.width * levelSet.mesh.height)

        # Record the time and area mismatch.
        times.append(runningTime)
        mismatches.append(mismatch)

        # Update the time of the next sample.
        nextSample += sampleInterval

        # Print statistics.
        print("%6.1f %11.4f" % (runningTime, mismatch))

        # Write level set and boundary segments to file.
        io.saveLevelSetVTK(len(times), levelSet)
        io.saveBoundarySegmentsTXT(len(times), boundary)

# Print results to file (distance vs time).
fileName = "shape-match_%5.4f.txt" % (temperature)
file = open(fileName, "w")
for i in range(0, len(times)):
    file.write("%lf %lf\n" % (times[i] - times[0], mismatches[i]))
file.close()

print("\nDone!")
