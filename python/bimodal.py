#  Copyright (c) 2015-2017 Lester Hedges <lester.hedges+pyslsm@gmail.com>
#                          Robert Jack   <r.jack@bath.ac.uk>
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
    bimodal.py

    An example showing perimeter minimisation with a shape matching
    constraint with a bimodal objective function.

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

    When a long trajectory is run the shape explores both minima within
    the dumbbell, with the time spent in each lobe proportional to the
    vertical offset between the lobe centres.

    The output file, "bimodal_*.txt", contains the measured x centre of mass,
    perimeter, mismatch, and centre of mass position vs time data for the
    optmisation run. Level set information for each sample interval is written
    to ParaView readable VTK files, "level-set_*.vtk". Boundary segment data
    is written to "boundary-segments_*.txt".
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
    if   y > uppergravityCutOff: return 1.0
    elif y < lowergravityCutOff: return reduce
    else:
        # Fractional distance along the gravity range.
        dy = (uppergravityCutOff - y) / gravityRange

        return (1.0 - dy*(1.0 - reduce));

# Boundary "centre of mass" function definition.
def computeCentreOfMass():
    # Zero variables.
    length = 0
    x = 0
    y = 0

    # Compute the "centre of mass" of the boundary.
    for i in range(0, len(boundary.points)):
        # Compute the weighted point length.
        ll = computePointLength(boundary.points[i])

        # Update variables.
        length += 0.5*ll;
        x += boundary.points[i].coord.x * 0.5*ll;
        y += boundary.points[i].coord.y * 0.5*ll;

    # Normalise.
    x /= length
    y /= length

    return x, y

# MAIN PROGRAM

# Boundary integral discretisation factor.
nDiscrete = 10

# Sensitivity reduction factor.
reduce = 0.65

# Maximum displacement per iteration, in units of the mesh spacing.
# This is the CFL limit.
moveLimit = 0.05

# Default temperature of the thermal bath.
temperature = 0.05

# Gravitational multiplier (must be between zero and 1).
gravityMult = 0.02

# Override temperature if command-line argument is passed.
if (len(sys.argv) > 1):
    temperature = float(sys.argv[1])

# Override sensitivity reduction factor if command-line argument is passed.
if (len(sys.argv) > 2):
    gravityMult = float(sys.argv[2])

# Set maximum running time.
maxTime = 1000000

# Set sampling interval.
sampleInterval = 1000

# Set time of the next sample (take an early snapshot).
nextSample = 1

# Width and height of the mesh.
width  = 70
height = 60;

# Hole vectors.
initialHoles = pyslsm.VectorHole()
targetHoles = pyslsm.VectorHole()

# Create a dumbbell from two horizontally offset holes.

# Centre the dumbell horizontally and vertically.
xCentre = 0.5*width
yCentre = 0.5*height

# Vertical offset between dumbbell lobes.
yOffset = 1

# Horizontal separation of lobe centres is 2.xOffset.
xOffset = 10

# Lobe radius.
radius = 20;

# Push holes into vectors.
targetHoles.append(pyslsm.Hole(xCentre+xOffset, yCentre+yOffset, radius))
targetHoles.append(pyslsm.Hole(xCentre-xOffset, yCentre, radius))

aa = xOffset / radius
area = (math.pi/2.0) + math.asin(aa) + aa*math.sqrt(1-aa*aa)
area *= 2.0 * radius * radius

# Set maximumum area mismatch (as fraction of total area).
minSizeOfShape = 0.2 / 3.0
minSizeOfShape *= width*height/area;

# Set the mismatch constraint.
maxMismatch = area/width/height * (1.0 - minSizeOfShape)

# Store dumbbell data (needed for sensitivity callback).
# They are used for implementing a cutoff for the "gravity".
uppergravityCutOff = targetHoles[0].coord.y + radius
lowergravityCutOff = targetHoles[0].coord.y - radius
gravityRange = uppergravityCutOff - lowergravityCutOff

# Physics stuff:
#    The "energy" of a boundary segment is l.g.y where l is length, g is gravity, y is height.
#    The original idea is that g.y=1      at uppergravityCutOff
#                            and  g.y=reduce at lowergravityCutOff
#    with linear interpolation in between
#        hence g.(upper-lower) = 1-reduce
#            g = (1-reduce)/(upper-lower) is the g-field
#        or
#            reduce = 1 - g.(upper-lower), which must be positive
#        hence
#            0 < g < 1/(upper-lower)

# The variable gravityMult is g.(upper-lower).
reduce = 1.0 - gravityMult

# Position and size of the initial trial circle.
initialHoleX = xCentre
initialHoleY = 0.5*yCentre
initialHoleRad = 0.5*radius

# Push the trial shape into the vector.
initialHoles.append(pyslsm.Hole(initialHoleX, initialHoleY, initialHoleRad))

# Initialise the level set domain.
levelSet = pyslsm.LevelSet(width, height, initialHoles, targetHoles, moveLimit, 6, True)

# Store the mesh area.
meshArea = levelSet.mesh.width * levelSet.mesh.height;

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
boundary.computeNormalVectors(levelSet);

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

# Centre of mass position measurements.
xPositions = pyslsm.VectorDouble()
yPositions = pyslsm.VectorDouble()

# Lambda values for the optimiser.
lambdas = pyslsm.VectorDouble([0, 0])

print("\nStarting bimodal demo...\n")

# Print output header.
print("----------------------------------------------------")
print("%9s %9s %10s %10s %10s" % ("Time", "Length", "Mismatch", "xPosition", "yPosition"))
print("----------------------------------------------------")

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

        # Compute the centre of mass of the shape.
        xPos, yPos = computeCentreOfMass()

        # Record the time, length, mismatch area, and positions.
        times.append(runningTime)
        lengths.append(length)
        mismatches.append(mismatch)
        xPositions.append(xPos-(width/2.0))
        yPositions.append(yPos-(width/2.0))

        # Update the time of the next sample.
        nextSample += sampleInterval

        # Print statistics.
        print("%9.4e %8.4f %10.4f %10.4f %10.4f" \
            % (runningTime, length, mismatch / meshArea, xPos-(width/2.0), yPos-(height/2.0)))

        # Write level set and boundary segments to file.
        io.saveLevelSetVTK(len(times), levelSet)
        io.saveBoundarySegmentsTXT(len(times), boundary)

# Print results to file.
fileName = "bimodal_%5.4f.txt" % (temperature)
file = open(fileName, "w")
for i in range(0, len(times)):
    file.write("%lf %lf %lf %lf %lf\n" \
        % (times[i] - times[0], lengths[i], mismatches[i] / meshArea, xPositions[i], yPositions[i]))
file.close()

print("\nDone!")
