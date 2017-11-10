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

""" minimise_area.py

    About:

    An example code showing a hole shrinking during unconstrained
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
"""

import math
import pyslsm

# Maximum displacement per iteration, in units of the mesh spacing.
# This is the CFL limit.
moveLimit = 0.5

# Set maximum running time.
maxTime = 50

# Set sampling interval.
sampleInterval = 1

# Set time of the next sample.
nextSample = 1

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

# Number of cycles since signed distance reinitialisation.
nReinit = 0

# Running time.
runningTime = 0

# Time measurements.
times = pyslsm.VectorDouble()

# Boundary length measurements.
lengths = pyslsm.VectorDouble()

# Lambda values for the optimiser.
lambdas = pyslsm.VectorDouble([0])

print("\nStarting unconstrained area minimisation demo...\n")

# Print output header.
print("---------------")
print("%6s %8s" % ("Time", "Length"))
print("---------------")

# Integrate until we exceed the maximum time.
while runningTime < maxTime:

    # Assign boundary point sensitivities.
    for i in range(0, len(boundary.points)):
        boundary.points[i].sensitivities[0] = 1.0

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
        # Record the time and boundary length.
        times.append(runningTime)
        lengths.append(boundary.length)

        # Update the time of the next sample.
        nextSample += sampleInterval

        # Print statistics.
        print("%6.1f %8.1f" % (runningTime, boundary.length))

        # Write level set and boundary segments to file.
        io.saveLevelSetVTK(len(times), levelSet)
        io.saveBoundarySegmentsTXT(len(times), boundary)

# Print results to file (distance vs time).
file = open("minimise_area.txt", "w")
for i in range(0, len(times)):
    file.write("%lf %lf\n" % (times[i] - times[0], (lengths[0] - lengths[i]) / (2 * math.pi)))
file.close()

print("\nDone!")
