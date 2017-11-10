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

"""
    dumbbell.py

    An example showing perimeter minimisation with a shape matching
    constraint.

    Here we construct a simple example system with two minima separated by
    a free energy barrier. The matched shape is a narrow-necked dumbbell
    constructed from two vertically offset, overlapping circles. The initial
    configuration is a circle centred in the upper lobe of the dumbbell.
    We construct two minima by making the perimeter (objective) sensitivities
    a function of the y coordinate of each boundary point. Sensitivities
    are linearly reduced between the centres of the upper and lower lobes,
    hence within the lower lobe it's possible to form a circle with a smaller
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
"""

import math
import pyslsm
import sys

# HELPER FUNCTIONS

# Constraint sensitivity function.
def computeConstraintSensitivity(coord, levelSet):
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

# Mismatch function.
def computeMismatch(mesh, targetArea):
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

# Perimeter function.
def computePerimeter(points):
    length = 0

    # Compute the total weighted boundary length.
    for i in range(0, len(points)):
        length += 0.5*computePointLength(points[i])

    return length

# Boundary point length function.
def computePointLength(point):
    length = 0

    # Store coordinates of the point.
    x1 = point.coord.x
    y1 = point.coord.y

    # Sum the distance to each neighbour.
    for i in range(0, point.nNeighbours):
        # Store coordinates of the neighbouring point.
        x2 = boundary.points[point.neighbours[i]].coord.x
        y2 = boundary.points[point.neighbours[i]].coord.y

        # Compute separation components.
        dx = x2 - x1
        dy = y2 - y1

        # Compute segment length.
        len = math.sqrt(dx*dx + dy*dy) / nDiscrete

        # Perform discrete boundary (line) integral.
        for j in range(0, nDiscrete):
            # Compute y position along segment.
            y = y1 + (j+0.5)*dy/nDiscrete

            # Add weighted segment length.
            length += computePerimeterWeight(y)*len

    return length

# Perimeter weight function.
def computePerimeterWeight(y):
    if   y > upperLobeCentre: return 1.0
    elif y < lowerLobeCentre: return reduce
    else:
        # Fractional distance from the centre of the upper dumbbell lobe.
        dy = (upperLobeCentre - y) / lobeSeparation

        return (1.0 - dy*(1.0 - reduce))

# MAIN PROGRAM

# Boundary integral discretisation factor.
nDiscrete = 10

# Sensitivity reduction factor.
reduce = 0.65

# Maximum displacement per iteration, in units of the mesh spacing.
# This is the CFL limit.
moveLimit = 0.05

# Default temperature of the thermal bath.
temperature = 0.02

# Override temperature if command-line argument is passed.
if (len(sys.argv) > 1):
    temperature = float(sys.argv[1])

# Override sensitivity reduction factor if command-line argument is passed.
if (len(sys.argv) > 2):
    reduce = float(sys.argv[2])

# Set maximum running time.
maxTime = 8000

# Set maximumum area mismatch.
maxMismatch = 0.2

# Set sampling interval.
sampleInterval = 40

# Set time of the next sample.
nextSample = 40

# Hole vectors.
initialHoles = pyslsm.VectorHole()
targetHoles = pyslsm.VectorHole()

# Create a dumbbell from two vertically offset holes.
targetHoles.append(pyslsm.Hole(50, 69, 20))
targetHoles.append(pyslsm.Hole(50, 31, 20))

# Store dumbbell data (needed for sensitivity callback).
upperLobeCentre = targetHoles[0].coord.y
lowerLobeCentre = targetHoles[1].coord.y
lobeSeparation = upperLobeCentre - lowerLobeCentre

# Initialise the system with a circle in the upper lobe.
initialHoles.append(pyslsm.Hole(50, 69, 15))

# Initialise a 100x100 level set domain.
levelSet = pyslsm.LevelSet(100, 100, initialHoles, targetHoles, moveLimit, 6, True)

# Store the mesh area.
meshArea = levelSet.mesh.width * levelSet.mesh.height

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
io.saveBoundarySegmentsTXT(0, boundary)

# Compute the element area fractions.
levelSet.computeAreaFractions(boundary)

# Store the target area fractions.
for i in range(0, levelSet.mesh.nElements):
    targetArea.append(levelSet.mesh.elements[i].area)

# Perform initial boundary discretisation.
boundary.discretise(levelSet)

# Compute the element area fractions.
levelSet.computeAreaFractions(boundary)

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

# Area mismatch measurements (constraint).
mismatches = pyslsm.VectorDouble()

# Lambda values for the optimiser.
lambdas = pyslsm.VectorDouble([0, 0])

print("Starting dumbbell demo...\n")

# Print output header.
print("--------------------------")
print("%6s %8s %10s" % ("Time", "Length", "Mismatch"))
print("--------------------------")

# Integrate until we exceed the maximum time.
while runningTime < maxTime:

    # Initialise the sensitivity object.
    sensitivity = pyslsm.Sensitivity()

    # Initialise the sensitivity callback.
    cb = pyslsm.Callback()
    cb.callback = computePointLength

    # Assign boundary point sensitivities.
    for i in range(0, len(boundary.points)):
        boundary.points[i].sensitivities[0] = \
            sensitivity.computeSensitivity(boundary.points[i], cb.callback)
        boundary.points[i].sensitivities[1] = \
            computeConstraintSensitivity(boundary.points[i].coord, levelSet)

    # Apply deterministic Ito correction.
    sensitivity.itoCorrection(boundary, temperature)

    # Time step associated with the iteration.
    timeStep = pyslsm.MutableFloat()

    # Constraint distance vector.
    constraintDistances = pyslsm.VectorDouble()

    # Current area mismatch.
    mismatch = computeMismatch(levelSet.mesh, targetArea)

    # Push current distance from constraint violation into vector.
    constraintDistances.append(meshArea*maxMismatch - mismatch)

    # Initialise the optimisation object.
    optimise = pyslsm.Optimise(boundary.points, constraintDistances, \
        lambdas, timeStep, levelSet.moveLimit)

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

    # Compute the element area fractions.
    levelSet.computeAreaFractions(boundary)

    # Compute the boundary point normal vectors.
    boundary.computeNormalVectors(levelSet)

    # Increment the time.
    runningTime += timeStep.value

    # Check if the next sample time has been reached.
    while runningTime >= nextSample:
        # Compute the weighted boundary perimeter.
        length = computePerimeter(boundary.points)

        # Record the time, length, and mismatch area.
        times.append(runningTime)
        lengths.append(length)
        mismatches.append(mismatch)

        # Update the time of the next sample.
        nextSample += sampleInterval

        # Print statistics.
        print("%6.1f %8.1f %10.4f" % (runningTime, length, mismatch / meshArea))

        # Write level set and boundary segments to file.
        io.saveLevelSetVTK(len(times), levelSet)
        io.saveBoundarySegmentsTXT(len(times), boundary)

# Print results to file.
fileName = "dumbbell_%5.4f.txt" % (temperature)
file = open(fileName, "w")
for i in range(0, len(times)):
    file.write("%lf %lf %lf\n" % (times[i] - times[0], lengths[i], mismatches[i] / meshArea))
file.close()

print("\nDone!")
