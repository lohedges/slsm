/*
  Copyright (c) 2015-2016-2016 Lester Hedges <lester.hedges+lsm@gmail.com>

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

#include "Boundary.h"
#include "LevelSet.h"

LevelSet::LevelSet(Mesh& mesh_, unsigned int bandWidth_) :
    nNodes(mesh_.nNodes),
    mesh(mesh_),
    bandWidth(bandWidth_)
{
    errno = EINVAL;
    check(bandWidth > 2, "Width of the narrow band must be greater than 2.");

    // Resize data structures.
    signedDistance.resize(nNodes);
    velocity.resize(nNodes);
    gradient.resize(nNodes);
    narrowBand.resize(nNodes);

    // Resize mines vector, 20% of nNodes is a reasonable estimate.
    mines.resize(int(0.2*nNodes));

    // Generate a Swiss cheese structure.
    initialise();

    // Initialise the narrow band.
    initialiseNarrowBand();

    return;

error:
    exit(EXIT_FAILURE);
}

LevelSet::LevelSet(Mesh& mesh_, unsigned int bandWidth_, const std::vector<Hole>& holes) :
    nNodes(mesh_.nNodes),
    mesh(mesh_),
    bandWidth(bandWidth_)
{
    errno = EINVAL;
    check(bandWidth > 2, "Width of the narrow band must be greater than 2.");

    // Resize data structures.
    signedDistance.resize(nNodes);
    velocity.resize(nNodes);
    gradient.resize(nNodes);
    narrowBand.resize(nNodes);

    // Resize mines vector, 20% of nNodes is a reasonable estimate.
    mines.resize(int(0.2*nNodes));

    // Initialise level set function from hole array.
    initialise(holes);

    // Initialise the narrow band.
    initialiseNarrowBand();

    return;

error:
    exit(EXIT_FAILURE);
}

void LevelSet::update()
{
}

void LevelSet::computeVelocities(const std::vector<BoundaryPoint>& boundaryPoints)
{
    // Initialise velocity (map boundary points to boundary nodes).
    initialiseVelocities(boundaryPoints);

    // Initialise fast marching method object.
    FastMarchingMethod fmm(mesh, false);

    // Re-initialise the signed distance function.
    fmm.march(signedDistance, velocity);
}

void LevelSet::computeGradients()
{
}

void LevelSet::reinitialise()
{
    // Initialise fast marching method object.
    FastMarchingMethod fmm(mesh, false);

    // Re-initialise the signed distance function.
    fmm.march(signedDistance);
}

void LevelSet::initialise()
{
    // Generate a swiss cheese arrangement of holes.
    // The holes have a default radius of 5.

    // Number of holes in x and y directions.
    unsigned int nx = std::round((double) mesh.width / 30);
    unsigned int ny = std::round((double) mesh.height / 30);

    // Number of holes.
    unsigned int n1 = (nx * ny);                // outer grid
    unsigned int n2 = ((nx - 1) * (ny - 1));    // inner grid
    unsigned int nHoles = n1 + n2;

    // Initialise a vector of holes.
    std::vector<Hole> holes(nHoles);

    // Hole separations.
    double dx, dy;

    // Check that mesh is large enough.
    check(((nx > 2) && (ny > 2)), "Mesh is too small for Swiss cheese initialisation.");

    // Initialise hole separations.
    dx = ((double) mesh.width / (2 * nx));
    dy = ((double) mesh.height / (2 * ny));

    // Calculate hole coordinates. (outer grid)
    for (unsigned int i=0;i<n1;i++)
    {
        // Work out x and y indices for hole.
        unsigned x = i % nx;
        unsigned int y = int(i / nx);

        // Set hole coordinates and radius.
        holes[i].coord.x = dx + (2 * x * dx);
        holes[i].coord.y = dy + (2 * y * dy);
        holes[i].r = 5;
    }

    // Calculate hole coordinates. (inner grid)
    for (unsigned int i=0;i<n2;i++)
    {
        // Work out x and y indices for hole.
        unsigned x = i % (nx - 1);
        unsigned int y = int(i / (nx -1));

        // Set hole coordinates and radius.
        holes[i + n1].coord.x = 2 * (dx + (x * dx));
        holes[i + n1].coord.y = 2 * (dy + (y * dy));
        holes[i + n1].r = 5;
    }

    // Now pass the holes to the initialise method.
    initialise(holes);

    return;

error:
    exit(EXIT_FAILURE);
}

void LevelSet::initialise(const std::vector<Hole>& holes)
{
    // First initialise LSF based on domain boundary.
    closestDomainBoundary();

    // Now test signed distance against the surface of each hole.
    // Update signed distance function when distance to hole surface
    // is less than the current value. Since this is only done once, we
    // use the simplest implementation possible.

    // Loop over all nodes.
    for (unsigned int i=0;i<nNodes;i++)
    {
        // Loop over all holes.
        for (unsigned int j=0;j<holes.size();j++)
        {
            // Work out x and y distance of the node from the hole centre.
            double dx = holes[j].coord.x - mesh.nodes[i].coord.x;
            double dy = holes[j].coord.y - mesh.nodes[i].coord.y;

            // Work out distance (Pythag).
            double dist = sqrt(dx*dx + dy*dy);

            // Signed distance from the hole surface.
            dist -= holes[j].r;

            // If distance is less than current value, then update.
            if (dist < signedDistance[i])
                signedDistance[i] = dist;
        }
    }
}

void LevelSet::closestDomainBoundary()
{
    // Initial LSF is distance from closest domain boundary.
    for (unsigned int i=0;i<nNodes;i++)
    {
        // Closest edge in x.
        unsigned int minX = std::min(mesh.nodes[i].coord.x, mesh.width - mesh.nodes[i].coord.x);

        // Closest edge in y.
        unsigned int minY = std::min(mesh.nodes[i].coord.y, mesh.height - mesh.nodes[i].coord.y);

        // Signed distance is the minimum of minX and minY;
        signedDistance[i] = double(std::min(minX, minY));
    }
}

void LevelSet::initialiseNarrowBand()
{
    unsigned int mineWidth = bandWidth - 1;

    // Reset the number of nodes in the narrow band.
    nNarrowBand = 0;

    // Reset the number of mines.
    nMines = 0;

    // Loop over all nodes.
    for (unsigned int i=0;i<nNodes;i++)
    {
        // Absolute value of the signed distance function.
        double absoluteSignedDistance = std::abs(signedDistance[i]);

        // Node lies inside band.
        if (absoluteSignedDistance < bandWidth)
        {
            // Flag node as active.
            mesh.nodes[i].isActive = true;

            // Update narrow band array.
            narrowBand[nNarrowBand] = i;

            // Increment number of nodes.
            nNarrowBand++;

            // Node lines at edge of band.
            if (absoluteSignedDistance > mineWidth)
            {
                // Node is a mine.
                mesh.nodes[i].isMine = true;

                // Update mine array.
                mines[nMines] = i;

                // Increment mine count.
                nMines++;

                // TODO:
                // Check when number of mines exceeds array size!
            }
        }
        // Node is outside band.
        else
        {
            // Flag node as inactive.
            mesh.nodes[i].isActive = false;
        }
    }
}

void LevelSet::initialiseVelocities(const std::vector<BoundaryPoint>& boundaryPoints)
{
    // Map boundary point velocities to nodes of the level set domain
    // using inverse squared distance interpolation.

    // Loop over all nodes.
    for (unsigned int i=0;i<nNodes;i++)
    {
        // Number of neighbouring boundary points.
        unsigned int nPoints = mesh.nodes[i].nBoundaryPoints;

        // Node has a single neighbouring boundary point.
        if (nPoints == 1)
            velocity[i] = boundaryPoints[mesh.nodes[i].boundaryPoints[0]].velocity;

        // Node has a two neighbouring boundary points.
        else if (nPoints == 2)
        {
            // The indices of the two boundary points.
            unsigned int point1 = mesh.nodes[i].boundaryPoints[0];
            unsigned int point2 = mesh.nodes[i].boundaryPoints[1];

            // Distances from the first point (x and y components).
            double dx1 = std::abs(mesh.nodes[i].coord.x - boundaryPoints[point1].coord.x);
            double dy1 = std::abs(mesh.nodes[i].coord.y - boundaryPoints[point1].coord.y);

            // Distances from the second point (x and y components).
            double dx2 = std::abs(mesh.nodes[i].coord.x - boundaryPoints[point2].coord.x);
            double dy2 = std::abs(mesh.nodes[i].coord.y - boundaryPoints[point2].coord.y);

            // Squared distances.
            double rSqd1 = dx1*dx1 + dy1*dy1;
            double rSqd2 = dx2*dx2 + dy2*dy2;

            // Weighting factors.
            double weight1 = 0.0;
            double weight2 = 0.0;

            // If one point lies exactly on a node, then use only that velocity.
            if (rSqd1 < 1e-6) weight1 = 1.0;
            else if (rSqd2 < 1e-6) weight2 = 1.0;

            // Determine weights (inverse squared distance).
            else
            {
                weight1 = 1.0 / rSqd1;
                weight2 = 1.0 / rSqd2;
            }

            // Calculate total weight.
            double totalWeight = weight1 + weight2;

            // Store interpolated nodal velocity: v = (w1*v1 + w2*v2) / (w1 + w2)
            velocity[i] = ((weight1 * boundaryPoints[point1].velocity)
                        + (weight2 * boundaryPoints[point2].velocity)) / totalWeight;
        }

        // No neighbouring points.
        else velocity[i] = 0;
    }
}
