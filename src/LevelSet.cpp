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

/*! \file LevelSet.cpp
    \brief A class for the level set function.
 */

namespace lsm
{
    LevelSet::LevelSet(Mesh& mesh_, double moveLimit_, unsigned int bandWidth_) :
        nNodes(mesh_.nNodes),
        moveLimit(moveLimit_),
        mesh(mesh_),
        bandWidth(bandWidth_)
    {
        int size = 0.2*mesh.nNodes;

        errno = EINVAL;
        lsm_check(bandWidth > 2, "Width of the narrow band must be greater than 2.");

        errno = EINVAL;
        lsm_check(((moveLimit > 0) && (moveLimit < 1)), "Move limit must be between 0 and 1.");

        // Resize data structures.
        signedDistance.resize(nNodes);
        velocity.resize(nNodes);
        gradient.resize(nNodes);
        narrowBand.resize(nNodes);

        // Make sure that memory is sufficient (for small test systems).
        size = std::max(25, size);
        mines.resize(size);

        // Generate a Swiss cheese structure.
        initialise();

        // Initialise the narrow band.
        initialiseNarrowBand();

        return;

    error:
        exit(EXIT_FAILURE);
    }

    LevelSet::LevelSet(Mesh& mesh_, const std::vector<Hole>& holes, double moveLimit_, unsigned int bandWidth_) :
        nNodes(mesh_.nNodes),
        moveLimit(moveLimit_),
        mesh(mesh_),
        bandWidth(bandWidth_)
    {
        int size = 0.2*mesh.nNodes;

        errno = EINVAL;
        lsm_check(bandWidth > 2, "Width of the narrow band must be greater than 2.");

        errno = EINVAL;
        lsm_check(((moveLimit > 0) && (moveLimit < 1)), "Move limit must be between 0 and 1.");

        // Resize data structures.
        signedDistance.resize(nNodes);
        velocity.resize(nNodes);
        gradient.resize(nNodes);
        narrowBand.resize(nNodes);

        // Make sure that memory is sufficient (for small test systems).
        size = std::max(25, size);
        mines.resize(size);

        // Initialise level set function from hole array.
        initialise(holes);

        // Initialise the narrow band.
        initialiseNarrowBand();

        return;

    error:
        exit(EXIT_FAILURE);
    }

    bool LevelSet::update(double timeStep)
    {
        // Loop over all nodes in the narrow band.
        for (unsigned int i=0;i<nNarrowBand;i++)
        {
            unsigned int node = narrowBand[i];
            signedDistance[node] -= timeStep * gradient[node] * velocity[node];
        }

        // Check mine nodes.
        for (unsigned int i=0;i<nMines;i++)
        {
            // Boundary is within one grid spacing of the mine.
            if (std::abs(signedDistance[mines[i]]) < 1.0)
            {
                // Reinitialise the signed distance function.
                reinitialise();

                return true;
            }
        }

        return false;
    }

    void LevelSet::computeVelocities(const std::vector<BoundaryPoint>& boundaryPoints)
    {
        // Initialise velocity (map boundary points to boundary nodes).
        initialiseVelocities(boundaryPoints);

        // Initialise fast marching method object.
        FastMarchingMethod fmm(mesh, false);

        // Reinitialise the signed distance function.
        fmm.march(signedDistance, velocity);
    }

    void LevelSet::computeGradients(const double timeStep)
    {
        // Compute gradient of the signed distance function using upwind finite difference.
        // This function assumes that velocities have already been calculated.

        // Reset gradients.
        std::fill(gradient.begin(), gradient.end(), 0.0);

        // Loop over all nodes in the narrow band region.
        for (unsigned int i=0;i<nNarrowBand;i++)
        {
            // Compute the nodal gradient.
            unsigned int node = narrowBand[i];
            gradient[node] = computeGradients(node, timeStep);
        }
    }

    void LevelSet::reinitialise()
    {
        // Initialise fast marching method object.
        FastMarchingMethod fmm(mesh, false);

        // Reinitialise the signed distance function.
        fmm.march(signedDistance);

        // Reinitialise the narrow band.
        initialiseNarrowBand();
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
        lsm_check(((nx > 2) && (ny > 2)), "Mesh is too small for Swiss cheese initialisation.");

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

                // Store interpolated nodal velocity: v = (w1*v1 + w2*v2]) / (w1 + w2)
                velocity[i] = ((weight1 * boundaryPoints[point1].velocity)
                            + (weight2 * boundaryPoints[point2].velocity)) / totalWeight;
            }

            // No neighbouring points.
            else velocity[i] = 0;
        }
    }

    double LevelSet::computeGradients(const unsigned int node, const double timeStep)
    {
        // Nodal coordinates.
        unsigned int x = mesh.nodes[node].coord.x;
        unsigned int y = mesh.nodes[node].coord.y;

        // Nodal signed distance.
        double lsf = signedDistance[node];

        // Whether gradient has been computed.
        bool isGradient = false;

        // Zero the gradient.
        double grad = 0;

        // Node is on the left edge.
        if (x == 0)
        {
            // Node is at bottom left corner.
            if (y == 0)
            {
                // If signed distance at nodes to right and above is the same, then use
                // the diagonal node for computing the gradient.
                if ((std::abs(signedDistance[mesh.xyToIndex[x+1][y]] - lsf) < 1e-6) &&
                    (std::abs(signedDistance[mesh.xyToIndex[x][y+1]] - lsf) < 1e-6))
                {
                    // Calculate signed distance to diagonal node.
                    grad = std::abs(lsf - signedDistance[mesh.xyToIndex[x+1][y+1]]);
                    grad *= sqrt(2.0);
                    isGradient = true;
                }
            }

            // Node is at top left corner.
            else if (y == mesh.height)
            {
                // If signed distance at nodes to right and below is the same, then use
                // the diagonal node for computing the gradient.
                if ((std::abs(signedDistance[mesh.xyToIndex[x+1][y]] - lsf) < 1e-6) &&
                    (std::abs(signedDistance[mesh.xyToIndex[x][y-1]] - lsf) < 1e-6))
                {
                    // Calculate signed distance to diagonal node.
                    grad = std::abs(lsf - signedDistance[mesh.xyToIndex[x+1][y-1]]);
                    grad *= sqrt(2.0);
                    isGradient = true;
                }
            }
        }

        // Node is on the right edge.
        else if (x == mesh.width)
        {
            // Node is at bottom right corner.
            if (y == 0)
            {
                // If signed distance at nodes to left and above is the same, then use
                // the diagonal node for computing the gradient.
                if ((std::abs(signedDistance[mesh.xyToIndex[x-1][y]] - lsf) < 1e-6) &&
                    (std::abs(signedDistance[mesh.xyToIndex[x][y+1]] - lsf) < 1e-6))
                {
                    // Calculate signed distance to diagonal node.
                    grad = std::abs(lsf - signedDistance[mesh.xyToIndex[x-1][y+1]]);
                    grad *= sqrt(2.0);
                    isGradient = true;
                }
            }

            // Node is at top right corner.
            else if (y == mesh.height)
            {
                // If signed distance at nodes to left and below is the same, then use
                // the diagonal node for computing the gradient.
                if ((std::abs(signedDistance[mesh.xyToIndex[x-1][y]] - lsf) < 1e-6) &&
                    (std::abs(signedDistance[mesh.xyToIndex[x][y-1]] - lsf) < 1e-6))
                {
                    // Calculate signed distance to diagonal node.
                    grad = std::abs(lsf - signedDistance[mesh.xyToIndex[x-1][y-1]]);
                    grad *= sqrt(2.0);
                    isGradient = true;
                }
            }
        }

        // Gradient hasn't already been calculated.
        if (!isGradient)
        {
            // Stencil values for the WENO approximation.
            double v1, v2, v3, v4, v5;

            int sign = velocity[node] < 0 ? -1 : 1;

            // Derivatives to right.

            // Node on left-hand edge.
            if (x == 0)
            {
                v5 = signedDistance[mesh.xyToIndex[3][y]] - signedDistance[mesh.xyToIndex[2][y]];
                v4 = signedDistance[mesh.xyToIndex[2][y]] - signedDistance[mesh.xyToIndex[1][y]];
                v3 = signedDistance[mesh.xyToIndex[1][y]] - signedDistance[mesh.xyToIndex[0][y]];

                // Approximate derivatives outside of domain.
                v2 = v3;
                v1 = v3;
            }

            // One node to right of left-hand edge.
            else if (x == 1)
            {
                v5 = signedDistance[mesh.xyToIndex[4][y]] - signedDistance[mesh.xyToIndex[3][y]];
                v4 = signedDistance[mesh.xyToIndex[3][y]] - signedDistance[mesh.xyToIndex[2][y]];
                v3 = signedDistance[mesh.xyToIndex[2][y]] - signedDistance[mesh.xyToIndex[1][y]];
                v2 = signedDistance[mesh.xyToIndex[1][y]] - signedDistance[mesh.xyToIndex[0][y]];

                // Approximate derivatives outside of domain.
                v1 = v2;
            }

            // Node on right-hand edge.
            else if (x == mesh.width)
            {
                v1 = signedDistance[mesh.xyToIndex[x-1][y]] - signedDistance[mesh.xyToIndex[x-2][y]];
                v2 = signedDistance[mesh.xyToIndex[x][y]] - signedDistance[mesh.xyToIndex[x-1][y]];

                // Approximate derivatives outside of domain.
                v3 = v2;
                v4 = v2;
                v5 = v2;
            }

            // One node to left of right-hand edge.
            else if (x == (mesh.width - 1))
            {
                v1 = signedDistance[mesh.xyToIndex[x-1][y]] - signedDistance[mesh.xyToIndex[x-2][y]];
                v2 = signedDistance[mesh.xyToIndex[x][y]] - signedDistance[mesh.xyToIndex[x-1][y]];
                v3 = signedDistance[mesh.xyToIndex[x+1][y]] - signedDistance[mesh.xyToIndex[x][y]];

                // Approximate derivatives outside of domain.
                v4 = v3;
                v5 = v3;
            }

            // Two nodes to left of right-hand edge.
            else if (x == (mesh.width - 2))
            {
                v1 = signedDistance[mesh.xyToIndex[x-1][y]] - signedDistance[mesh.xyToIndex[x-2][y]];
                v2 = signedDistance[mesh.xyToIndex[x][y]] - signedDistance[mesh.xyToIndex[x-1][y]];
                v3 = signedDistance[mesh.xyToIndex[x+1][y]] - signedDistance[mesh.xyToIndex[x][y]];
                v4 = signedDistance[mesh.xyToIndex[x+2][y]] - signedDistance[mesh.xyToIndex[x+1][y]];

                // Approximate derivatives outside of domain.
                v5 = v4;
            }

            // Node lies in bulk.
            else
            {
                v1 = signedDistance[mesh.xyToIndex[x-1][y]] - signedDistance[mesh.xyToIndex[x-2][y]];
                v2 = signedDistance[mesh.xyToIndex[x][y]] - signedDistance[mesh.xyToIndex[x-1][y]];
                v3 = signedDistance[mesh.xyToIndex[x+1][y]] - signedDistance[mesh.xyToIndex[x][y]];
                v4 = signedDistance[mesh.xyToIndex[x+2][y]] - signedDistance[mesh.xyToIndex[x+1][y]];
                v5 = signedDistance[mesh.xyToIndex[x+3][y]] - signedDistance[mesh.xyToIndex[x+2][y]];
            }

            double gradRight = sign * gradWENO(v1, v2, v3, v4, v5);

            // Derivatives to left.

            // Node on right-hand edge.
            if (x == mesh.width)
            {
                v5 = signedDistance[mesh.xyToIndex[x-2][y]] - signedDistance[mesh.xyToIndex[x-3][y]];
                v4 = signedDistance[mesh.xyToIndex[x-1][y]] - signedDistance[mesh.xyToIndex[x-2][y]];
                v3 = signedDistance[mesh.xyToIndex[x][y]] - signedDistance[mesh.xyToIndex[x-1][y]];

                // Approximate derivatives outside of domain.
                v2 = v3;
                v1 = v3;
            }

            // One node to left of right-hand edge.
            else if (x == (mesh.width-1))
            {
                v5 = signedDistance[mesh.xyToIndex[x-2][y]] - signedDistance[mesh.xyToIndex[x-3][y]];
                v4 = signedDistance[mesh.xyToIndex[x-1][y]] - signedDistance[mesh.xyToIndex[x-2][y]];
                v3 = signedDistance[mesh.xyToIndex[x][y]] - signedDistance[mesh.xyToIndex[x-1][y]];
                v2 = signedDistance[mesh.xyToIndex[x+1][y]] - signedDistance[mesh.xyToIndex[x][y]];

                // Approximate derivatives outside of domain.
                v1 = v2;
            }

            // Node on left-hand edge.
            else if (x == 0)
            {
                v1 = signedDistance[mesh.xyToIndex[2][y]] - signedDistance[mesh.xyToIndex[1][y]];
                v2 = signedDistance[mesh.xyToIndex[1][y]] - signedDistance[mesh.xyToIndex[0][y]];

                // Approximate derivatives outside of domain.
                v3 = v2;
                v4 = v2;
                v5 = v2;
            }

            // One node to right of left-hand edge.
            else if (x == 1)
            {
                v1 = signedDistance[mesh.xyToIndex[3][y]] - signedDistance[mesh.xyToIndex[2][y]];
                v2 = signedDistance[mesh.xyToIndex[2][y]] - signedDistance[mesh.xyToIndex[1][y]];
                v3 = signedDistance[mesh.xyToIndex[1][y]] - signedDistance[mesh.xyToIndex[0][y]];

                // Approximate derivatives outside of domain.
                v4 = v3;
                v5 = v3;
            }

            // Two nodes to right of left-hand edge.
            else if (x == 2)
            {
                v1 = signedDistance[mesh.xyToIndex[4][y]] - signedDistance[mesh.xyToIndex[3][y]];
                v2 = signedDistance[mesh.xyToIndex[3][y]] - signedDistance[mesh.xyToIndex[2][y]];
                v3 = signedDistance[mesh.xyToIndex[2][y]] - signedDistance[mesh.xyToIndex[1][y]];
                v4 = signedDistance[mesh.xyToIndex[1][y]] - signedDistance[mesh.xyToIndex[0][y]];

                // Approximate derivatives outside of domain.
                v5 = v4;
            }

            // Node lies in bulk.
            else
            {
                v1 = signedDistance[mesh.xyToIndex[x+2][y]] - signedDistance[mesh.xyToIndex[x+1][y]];
                v2 = signedDistance[mesh.xyToIndex[x+1][y]] - signedDistance[mesh.xyToIndex[x][y]];
                v3 = signedDistance[mesh.xyToIndex[x][y]] - signedDistance[mesh.xyToIndex[x-1][y]];
                v4 = signedDistance[mesh.xyToIndex[x-1][y]] - signedDistance[mesh.xyToIndex[x-2][y]];
                v5 = signedDistance[mesh.xyToIndex[x-2][y]] - signedDistance[mesh.xyToIndex[x-3][y]];
            }

            double gradLeft = sign * gradWENO(v1, v2, v3, v4, v5);

            // Upward derivatives.

            // Node on bottom edge.
            if (y == 0)
            {
                v5 = signedDistance[mesh.xyToIndex[x][3]] - signedDistance[mesh.xyToIndex[x][2]];
                v4 = signedDistance[mesh.xyToIndex[x][2]] - signedDistance[mesh.xyToIndex[x][1]];
                v3 = signedDistance[mesh.xyToIndex[x][1]] - signedDistance[mesh.xyToIndex[x][0]];

                // Approximate derivatives outside of domain.
                v2 = v3;
                v1 = v3;
            }

            // One node above bottom edge.
            else if (y == 1)
            {
                v5 = signedDistance[mesh.xyToIndex[x][4]] - signedDistance[mesh.xyToIndex[x][3]];
                v4 = signedDistance[mesh.xyToIndex[x][3]] - signedDistance[mesh.xyToIndex[x][2]];
                v3 = signedDistance[mesh.xyToIndex[x][2]] - signedDistance[mesh.xyToIndex[x][1]];
                v2 = signedDistance[mesh.xyToIndex[x][1]] - signedDistance[mesh.xyToIndex[x][0]];

                // Approximate derivatives outside of domain.
                v1 = v2;
            }

            // Node is on top edge.
            else if (y == mesh.height)
            {
                v1 = signedDistance[mesh.xyToIndex[x][y-1]] - signedDistance[mesh.xyToIndex[x][y-2]];
                v2 = signedDistance[mesh.xyToIndex[x][y]] - signedDistance[mesh.xyToIndex[x][y-1]];

                // Approximate derivatives outside of domain.
                v3 = v2;
                v4 = v2;
                v5 = v2;
            }

            // One node below top edge.
            else if (y == (mesh.height - 1))
            {
                v1 = signedDistance[mesh.xyToIndex[x][y-1]] - signedDistance[mesh.xyToIndex[x][y-2]];
                v2 = signedDistance[mesh.xyToIndex[x][y]] - signedDistance[mesh.xyToIndex[x][y-1]];
                v3 = signedDistance[mesh.xyToIndex[x][y+1]] - signedDistance[mesh.xyToIndex[x][y]];

                // Approximate derivatives outside of domain.
                v4 = v3;
                v5 = v3;
            }

            // Two nodes below top edge.
            else if (y == (mesh.height - 2))
            {
                v1 = signedDistance[mesh.xyToIndex[x][y-1]] - signedDistance[mesh.xyToIndex[x][y-2]];
                v2 = signedDistance[mesh.xyToIndex[x][y]] - signedDistance[mesh.xyToIndex[x][y-1]];
                v3 = signedDistance[mesh.xyToIndex[x][y+1]] - signedDistance[mesh.xyToIndex[x][y]];
                v4 = signedDistance[mesh.xyToIndex[x][y+2]] - signedDistance[mesh.xyToIndex[x][y+1]];

                // Approximate derivatives outside of domain.
                v5 = v4;
            }

            // Node lies in bulk.
            else
            {
                v1 = signedDistance[mesh.xyToIndex[x][y-1]] - signedDistance[mesh.xyToIndex[x][y-2]];
                v2 = signedDistance[mesh.xyToIndex[x][y]] - signedDistance[mesh.xyToIndex[x][y-1]];
                v3 = signedDistance[mesh.xyToIndex[x][y+1]] - signedDistance[mesh.xyToIndex[x][y]];
                v4 = signedDistance[mesh.xyToIndex[x][y+2]] - signedDistance[mesh.xyToIndex[x][y+1]];
                v5 = signedDistance[mesh.xyToIndex[x][y+3]] - signedDistance[mesh.xyToIndex[x][y+2]];
            }

            double gradUp = sign * gradWENO(v1, v2, v3, v4, v5);

            // Downward derivative.

            // Node on right edge.
            if (y == mesh.width)
            {
                v5 = signedDistance[mesh.xyToIndex[x][y-2]] - signedDistance[mesh.xyToIndex[x][y-3]];
                v4 = signedDistance[mesh.xyToIndex[x][y-1]] - signedDistance[mesh.xyToIndex[x][y-2]];
                v3 = signedDistance[mesh.xyToIndex[x][y]] - signedDistance[mesh.xyToIndex[x][y-1]];

                // Approximate derivatives outside of domain.
                v2 = v3;
                v1 = v3;
            }

            // One node left of right edge.
            else if (y == (mesh.width - 1))
            {
                v5 = signedDistance[mesh.xyToIndex[x][y-2]] - signedDistance[mesh.xyToIndex[x][y-3]];
                v4 = signedDistance[mesh.xyToIndex[x][y-1]] - signedDistance[mesh.xyToIndex[x][y-2]];
                v3 = signedDistance[mesh.xyToIndex[x][y]] - signedDistance[mesh.xyToIndex[x][y-1]];
                v2 = signedDistance[mesh.xyToIndex[x][y+1]] - signedDistance[mesh.xyToIndex[x][y]];

                // Approximate derivatives outside of domain.
                v1 = v2;
            }

            // Node lies on bottom edge
            else if (y == 0)
            {
                v1 = signedDistance[mesh.xyToIndex[x][2]] - signedDistance[mesh.xyToIndex[x][1]];
                v2 = signedDistance[mesh.xyToIndex[x][1]] - signedDistance[mesh.xyToIndex[x][0]];

                // Approximate derivatives outside of domain.
                v3 = v2;
                v4 = v2;
                v5 = v2;
            }

            // One node above bottom edge.
            else if (y == 1)
            {
                v1 = signedDistance[mesh.xyToIndex[x][3]] - signedDistance[mesh.xyToIndex[x][2]];
                v2 = signedDistance[mesh.xyToIndex[x][2]] - signedDistance[mesh.xyToIndex[x][1]];
                v3 = signedDistance[mesh.xyToIndex[x][1]] - signedDistance[mesh.xyToIndex[x][0]];

                // Approximate derivatives outside of domain.
                v4 = v3;
                v5 = v3;
            }

            // Two nodes above bottom edge.
            else if (y == 2)
            {
                v1 = signedDistance[mesh.xyToIndex[x][4]] - signedDistance[mesh.xyToIndex[x][3]];
                v2 = signedDistance[mesh.xyToIndex[x][3]] - signedDistance[mesh.xyToIndex[x][2]];
                v3 = signedDistance[mesh.xyToIndex[x][2]] - signedDistance[mesh.xyToIndex[x][1]];
                v4 = signedDistance[mesh.xyToIndex[x][1]] - signedDistance[mesh.xyToIndex[x][0]];

                // Approximate derivatives outside of domain.
                v5 = v4;
            }

            // Node lies in bulk.
            else
            {
                v1 = signedDistance[mesh.xyToIndex[x][y-2]] - signedDistance[mesh.xyToIndex[x][y-3]];
                v2 = signedDistance[mesh.xyToIndex[x][y-1]] - signedDistance[mesh.xyToIndex[x][y-2]];
                v3 = signedDistance[mesh.xyToIndex[x][y]] - signedDistance[mesh.xyToIndex[x][y-1]];
                v4 = signedDistance[mesh.xyToIndex[x][y+1]] - signedDistance[mesh.xyToIndex[x][y]];
                v5 = signedDistance[mesh.xyToIndex[x][y+2]] - signedDistance[mesh.xyToIndex[x][y+1]];
            }

            double gradDown = sign * gradWENO(v1, v2, v3, v4, v5);

            // Compute gradient using upwind scheme.

            if (gradDown > 0)   grad += gradDown * gradDown;
            if (gradLeft > 0)   grad += gradLeft * gradLeft;
            if (gradUp < 0)     grad += gradUp * gradUp;
            if (gradRight < 0)  grad += gradRight * gradRight;

            grad = sqrt(grad);
        }

        // Check that gradient doesn't take the boundary outside of the domain.
        if ((x == 0) || (x == mesh.width) || (y == 0) || (y == mesh.height))
        {
            // Updated signed distance is positive.
            if ((lsf - (timeStep * grad * velocity[node])) > 0)
                grad = lsf / (timeStep * velocity[node]);
        }

        // Return gradient.
        return grad;
    }

    double LevelSet::gradWENO(double v1, double v2, double v3, double v4, double v5)
    {
        // Approximate the gradient using the 5th order WENO approximation.
        // See: http://www.scholarpedia.org/article/WENO_methods

        double oneThird     = 1.0 / 3.0;
        double oneQuarter   = 1.0 / 4.0;
        double oneEigth     = 1.0 / 8.0;
        double oneSixteenth = 1.0 / 16.0;

        // Compute the beta values for each stencil.
        double beta1 = oneThird * (4*v1*v1 - 19*v1*v2 + 25*v2*v2 + 11*v1*v3 - 31*v2*v3 + 10*v3*v3);
        double beta2 = oneThird * (4*v2*v2 - 13*v2*v3 + 13*v3*v3 + 5*v2*v4 - 13*v3*v4 + 4*v4*v4);
        double beta3 = oneThird * (10*v3*v3 - 31*v3*v4 + 25*v4*v4 + 11*v3*v5 - 19*v4*v5 + 4*v5*v5);

        // Complete the non-linear weights.
        double w1 = oneSixteenth / ((1e-6 + beta1) * (1e-6 + beta1));
        double w2 = (5 * oneEigth) / ((1e-6 + beta2) * (1e-6 + beta2));
        double w3 = (5 * oneSixteenth) / ((1e-6 + beta3) * (1e-6 + beta3));

        // Normalise weights.
        double totalWeight = w1 + w2 + w3;
        w1 /= totalWeight;
        w2 /= totalWeight;
        w3 /= totalWeight;

        // Sum three stencil components.
        double grad = w1 * ((3 * oneEigth * v1) - (5 * oneQuarter * v2) + (15 * oneEigth * v3))
                    + w2 * ((-oneEigth * v2) + (3 * oneQuarter * v3) + (3 * oneEigth * v4))
                    + w3 * ((3 * oneEigth * v3) + (3 * oneQuarter * v4) - (oneEigth * v5));

        return grad;
    }
}
