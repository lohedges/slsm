/*
  Copyright (c) 2015-2016 Lester Hedges <lester.hedges+slsm@gmail.com>

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

#include <cmath>

#include "Boundary.h"
#include "LevelSet.h"
#include "Mesh.h"

/*! \file Boundary.cpp
    \brief A class for the discretised boundary.
 */

namespace slsm
{
    Boundary::Boundary()
    {
    }

    void Boundary::discretise(LevelSet& levelSet, bool isTarget)
    {
        // Clear and reserve vector memory.
        points.clear();
        segments.clear();
        points.reserve(levelSet.mesh.nNodes);
        segments.reserve(levelSet.mesh.nNodes);

        // Reset the number of points and segments.
        nPoints = nSegments = 0;

        // Zero the boundary length.
        length = 0;

        // Initialise a pointer to the signed distance vector.
        std::vector<double>* signedDistance;

        // Point to the target signed distance.
        if (isTarget) signedDistance = &levelSet.target;

        // Point to the current signed distance function.
        else signedDistance = &levelSet.signedDistance;

        // Compute the status of nodes and elements in level-set mesh.
        computeMeshStatus(levelSet.mesh, signedDistance);

        // Loop over all elements.
        for (unsigned int i=0;i<levelSet.mesh.nElements;i++)
        {
            // The element isn't outside of the structure.
            if (levelSet.mesh.elements[i].status != ElementStatus::OUTSIDE)
            {
                // Number of cut edges.
                unsigned int nCut = 0;

                // Existing boundary points associated with a node.
                unsigned int boundaryPoints[4];

                // Look at each edge of the node to determine whether it is cut,
                // or if it's part of the boundary.
                for (unsigned int j=0;j<4;j++)
                {
                    // Index of first node on edge.
                    unsigned int n1 = levelSet.mesh.elements[i].nodes[j];

                    // Index of second node (reconnecting to 0th node).
                    // Edge connectivity goes: 0 --> 1, 1 --> 2, 2 --> 3, 3 --> 0
                    unsigned int n2 = (j == 3) ? 0 : (j + 1);

                    // Convert to node index.
                    n2 = levelSet.mesh.elements[i].nodes[n2];

                    // If not performing discretisation of a target structure check that at least
                    // one node lies in the narrow band region, or is masked.
                    if (isTarget ||
                    (levelSet.mesh.nodes[n1].isActive | levelSet.mesh.nodes[n1].isMasked) ||
                    (levelSet.mesh.nodes[n2].isActive | levelSet.mesh.nodes[n2].isMasked))
                    {
                        // One node is inside, the other is outside. The edge is cut.
                        if ((levelSet.mesh.nodes[n1].status|levelSet.mesh.nodes[n2].status) == NodeStatus::CUT)
                        {
                            // Compute the distance from node 1 to the intersection point (by interpolation).
                            double d = (*signedDistance)[n1]
                                     / ((*signedDistance)[n1] - (*signedDistance)[n2]);

                            // Initialise boundary point coordinate.
                            Coord coord;

                            // Make sure that the boundary point hasn't already been added.
                            int index = isAdded(levelSet.mesh, coord, n1, j, d);

                            // Boundary point is new.
                            if (index < 0)
                            {
                                levelSet.mesh.nodes[n1].boundaryPoints[levelSet.mesh.nodes[n1].nBoundaryPoints] = nPoints;
                                levelSet.mesh.nodes[n2].boundaryPoints[levelSet.mesh.nodes[n2].nBoundaryPoints] = nPoints;
                                levelSet.mesh.nodes[n1].nBoundaryPoints++;
                                levelSet.mesh.nodes[n2].nBoundaryPoints++;

                                // Store boundary point for cut edge.
                                boundaryPoints[nCut] = nPoints;

                                // Initialise boundary point.
                                BoundaryPoint point;
                                initialisePoint(levelSet, point, coord);
                                points.push_back(point);

                                // Increment number of boundary points.
                                nPoints++;
                            }
                            else
                            {
                                // Store existing boundary point.
                                boundaryPoints[nCut] = index;
                            }

                            // Increment number of cut edges.
                            nCut++;
                        }

                        // Both nodes lie on the boundary.
                        else if ((levelSet.mesh.nodes[n1].status & NodeStatus::BOUNDARY) &&
                            (levelSet.mesh.nodes[n2].status & NodeStatus::BOUNDARY))
                        {
                            // Initialise boundary point coordinate.
                            Coord coord;

                            // Create boundary segment.
                            BoundarySegment segment;

                            // Assign element index.
                            segment.element = i;

                            // Make sure that the start boundary point hasn't already been added.
                            int index = isAdded(levelSet.mesh, coord, n1, 0, 0);

                            // Boundary point is new.
                            if (index < 0)
                            {
                                // Set index equal to current number of points.
                                index = nPoints;

                                levelSet.mesh.nodes[n1].boundaryPoints[levelSet.mesh.nodes[n1].nBoundaryPoints] = nPoints;
                                levelSet.mesh.nodes[n1].nBoundaryPoints++;

                                // Initialise boundary point.
                                BoundaryPoint point;
                                initialisePoint(levelSet, point, coord);
                                points.push_back(point);

                                // Increment number of boundary points.
                                nPoints++;
                            }

                            // Assign start point index.
                            segment.start = index;

                            // Make sure that the end boundary point hasn't already been added.
                            index = isAdded(levelSet.mesh, coord, n2, 0, 0);

                            // Boundary point is new.
                            if (index < 0)
                            {
                                // Set index equal to current number of points.
                                index = nPoints;

                                levelSet.mesh.nodes[n2].boundaryPoints[levelSet.mesh.nodes[n2].nBoundaryPoints] = nPoints;
                                levelSet.mesh.nodes[n2].nBoundaryPoints++;

                                // Initialise boundary point.
                                BoundaryPoint point;
                                initialisePoint(levelSet, point, coord);
                                points.push_back(point);

                                // Increment number of boundary points.
                                nPoints++;
                            }

                            // Assign end point index.
                            segment.end = index;

                            // Compute the length of the boundary segment.
                            segment.length = segmentLength(segment);

                            // Update total boundary length.
                            length += segment.length;

                            // Create element to segment lookup.
                            levelSet.mesh.elements[i].boundarySegments[levelSet.mesh.elements[i].nBoundarySegments] = nSegments;
                            levelSet.mesh.elements[i].nBoundarySegments++;

                            // Add segment to vector.
                            segments.push_back(segment);
                            nSegments++;
                        }
                    }
                }

                // If the element was cut, determine the boundary segment(s).

                // If two edges are cut, then a boundary segment must cross both.
                if (nCut == 2)
                {
                    // Create boundary segment.
                    BoundarySegment segment;
                    segment.start = boundaryPoints[0];
                    segment.end = boundaryPoints[1];
                    segment.element = i;

                    // Compute the length of the boundary segment.
                    segment.length = segmentLength(segment);

                    // Update total boundary length.
                    length += segment.length;

                    // Create element to segment lookup.
                    levelSet.mesh.elements[i].boundarySegments[levelSet.mesh.elements[i].nBoundarySegments] = nSegments;
                    levelSet.mesh.elements[i].nBoundarySegments++;

                    // Add segment to vector.
                    segments.push_back(segment);
                    nSegments++;
                }

                // If there is only one cut edge, then the boundary must also cross an element node.
                else if (nCut == 1)
                {
                    // Find a node that is on the boundary and has a neighbour that is outside.
                    for (unsigned int j=0;j<4;j++)
                    {
                        // Node index.
                        unsigned int node = levelSet.mesh.elements[i].nodes[j];

                        // Node is on the boundary, check its neighbours.
                        if (levelSet.mesh.nodes[node].status & NodeStatus::BOUNDARY)
                        {
                            // Index of next node.
                            unsigned int nAfter = (j == 3) ? 0 : (j + 1);

                            // Index of previous node.
                            unsigned int nBefore = (j == 0) ? 3 : (j - 1);

                            // Convert to node indices.
                            nAfter = levelSet.mesh.elements[i].nodes[nAfter];
                            nBefore = levelSet.mesh.elements[i].nodes[nBefore];

                            // If a neighbour is outside the boundary, then add a boundary segment.
                            if ((levelSet.mesh.nodes[nAfter].status & NodeStatus::OUTSIDE) ||
                                (levelSet.mesh.nodes[nBefore].status & NodeStatus::OUTSIDE))
                            {
                                // Create boundary segment.
                                BoundarySegment segment;
                                segment.start = boundaryPoints[0];
                                segment.element = i;

                                // Initialise boundary point coordinate.
                                Coord coord;

                                // Make sure that the end boundary point hasn't already been added.
                                int index = isAdded(levelSet.mesh, coord, node, 0, 0);

                                // Boundary point is new.
                                if (index < 0)
                                {
                                    // Set index equal to current number of points.
                                    index = nPoints;

                                    levelSet.mesh.nodes[node].boundaryPoints[levelSet.mesh.nodes[node].nBoundaryPoints] = nPoints;
                                    levelSet.mesh.nodes[node].nBoundaryPoints++;

                                    // Initialise boundary point.
                                    BoundaryPoint point;
                                    initialisePoint(levelSet, point, coord);
                                    points.push_back(point);

                                    // Increment number of boundary points.
                                    nPoints++;
                                }

                                // Assign end point index.
                                segment.end = index;

                                // Compute the length of the boundary segment.
                                segment.length = segmentLength(segment);

                                // Update total boundary length.
                                length += segment.length;

                                // Create element to segment lookup.
                                levelSet.mesh.elements[i].boundarySegments[levelSet.mesh.elements[i].nBoundarySegments] = nSegments;
                                levelSet.mesh.elements[i].nBoundarySegments++;

                                // Add segment to vector.
                                segments.push_back(segment);
                                nSegments++;
                            }
                        }
                    }
                }

                // If there are four cut edges, then determine which
                // boundary node pairs form the boundary.
                else if (nCut == 4)
                {
                    double lsfSum = 0;

                    // Evaluate level set value at element centre.
                    for (unsigned int j=0;j<4;j++)
                    {
                        // Node index.
                        unsigned int node = levelSet.mesh.elements[i].nodes[j];

                        lsfSum += (*signedDistance)[node];
                    }

                    // Create boundary segment.
                    BoundarySegment segment;

                    // Store the status of the first node.
                    unsigned int node = levelSet.mesh.elements[i].nodes[0];
                    NodeStatus::NodeStatus status = levelSet.mesh.nodes[node].status;

                    if (((status & NodeStatus::INSIDE) && (lsfSum > 0)) ||
                        ((status & NodeStatus::OUTSIDE) && (lsfSum < 0)))
                    {
                        segment.start = boundaryPoints[0];
                        segment.end = boundaryPoints[1];
                        segment.element = i;

                        // Compute the length of the boundary segment.
                        segment.length = segmentLength(segment);

                        // Update total boundary length.
                        length += segment.length;

                        // Create element to segment lookup.
                        levelSet.mesh.elements[i].boundarySegments[levelSet.mesh.elements[i].nBoundarySegments] = nSegments;
                        levelSet.mesh.elements[i].nBoundarySegments++;

                        // Add segment to vector.
                        segments.push_back(segment);
                        nSegments++;

                        segment.start = boundaryPoints[2];
                        segment.end = boundaryPoints[3];
                        segment.element = i;

                        // Compute the length of the boundary segment.
                        segment.length = segmentLength(segment);

                        // Update total boundary length.
                        length += segment.length;

                        // Create element to segment lookup.
                        levelSet.mesh.elements[i].boundarySegments[levelSet.mesh.elements[i].nBoundarySegments] = nSegments;
                        levelSet.mesh.elements[i].nBoundarySegments++;

                        // Add segment to vector.
                        segments.push_back(segment);
                        nSegments++;
                    }

                    else
                    {
                        segment.start = boundaryPoints[0];
                        segment.end = boundaryPoints[3];
                        segment.element = i;

                        // Compute the length of the boundary segment.
                        segment.length = segmentLength(segment);

                        // Update total boundary length.
                        length += segment.length;

                        // Create element to segment lookup.
                        levelSet.mesh.elements[i].boundarySegments[levelSet.mesh.elements[i].nBoundarySegments] = nSegments;
                        levelSet.mesh.elements[i].nBoundarySegments++;

                        // Add segment to vector.
                        segments.push_back(segment);
                        nSegments++;

                        segment.start = boundaryPoints[1];
                        segment.end = boundaryPoints[2];
                        segment.element = i;

                        // Compute the length of the boundary segment.
                        segment.length = segmentLength(segment);

                        // Update total boundary length.
                        length += segment.length;

                        // Create element to segment lookup.
                        levelSet.mesh.elements[i].boundarySegments[levelSet.mesh.elements[i].nBoundarySegments] = nSegments;
                        levelSet.mesh.elements[i].nBoundarySegments++;

                        // Add segment to vector.
                        segments.push_back(segment);
                        nSegments++;
                    }

                    // Update element status to indicate whether centre is in or out.
                    levelSet.mesh.elements[i].status = (lsfSum > 0) ? ElementStatus::CENTRE_INSIDE : ElementStatus::CENTRE_OUTSIDE;
                }

                // If no edges are cut and element is not inside structure
                // then the boundary segment must cross the diagonal.
                else if ((nCut == 0) && (levelSet.mesh.elements[i].status != ElementStatus::INSIDE))
                {
                    // Node index.
                    unsigned int node;

                    // Find the two boundary nodes.
                    for (unsigned int j=0;j<4;j++)
                    {
                        // Node index.
                        node = levelSet.mesh.elements[i].nodes[j];

                        if (levelSet.mesh.nodes[node].status & NodeStatus::BOUNDARY)
                        {
                            boundaryPoints[nCut] = node;
                            nCut++;
                        }
                    }

                    // Create boundary segment.
                    BoundarySegment segment;
                    segment.element = i;

                    // Initialise boundary point coordinate.
                    Coord coord;

                    // Make sure that the start boundary point hasn't already been added.
                    node = boundaryPoints[0];
                    int index = isAdded(levelSet.mesh, coord, node, 0, 0);

                    // Boundary point is new.
                    if (index < 0)
                    {
                        // Set index equal to current number of points.
                        index = nPoints;

                        levelSet.mesh.nodes[node].boundaryPoints[levelSet.mesh.nodes[node].nBoundaryPoints] = nPoints;
                        levelSet.mesh.nodes[node].nBoundaryPoints++;

                        // Initialise boundary point.
                        BoundaryPoint point;
                        initialisePoint(levelSet, point, coord);
                        points.push_back(point);

                        // Increment number of boundary points.
                        nPoints++;
                    }

                    // Assign start point index.
                    segment.start = index;

                    // Make sure that the end boundary point hasn't already been added.
                    node = boundaryPoints[1];
                    index = isAdded(levelSet.mesh, coord, node, 0, 0);

                    // Boundary point is new.
                    if (index < 0)
                    {
                        // Set index equal to current number of points.
                        index = nPoints;

                        levelSet.mesh.nodes[node].boundaryPoints[levelSet.mesh.nodes[node].nBoundaryPoints] = nPoints;
                        levelSet.mesh.nodes[node].nBoundaryPoints++;

                        // Initialise boundary point.
                        BoundaryPoint point;
                        initialisePoint(levelSet, point, coord);
                        points.push_back(point);

                        // Increment number of boundary points.
                        nPoints++;
                    }

                    // Assign end point index.
                    segment.end = index;

                    // Compute the length of the boundary segment.
                    segment.length = segmentLength(segment);

                    // Update total boundary length.
                    length += segment.length;

                    // Create element to segment lookup.
                    levelSet.mesh.elements[i].boundarySegments[levelSet.mesh.elements[i].nBoundarySegments] = nSegments;
                    levelSet.mesh.elements[i].nBoundarySegments++;

                    // Add segment to vector.
                    segments.push_back(segment);
                    nSegments++;
                }
            }
        }

        // Work out boundary integral length associated with each boundary point.
        computePointLengths();
    }

    void Boundary::computeNormalVectors(const LevelSet& levelSet)
    {
        // Whether the normal vector at a boundary point has been set.
        bool isSet[nPoints];

        // Weighting factor for each point.
        double weight[nPoints];

        // Initialise arrays.
        for (unsigned int i=0;i<nPoints;i++)
        {
            isSet[i] = false;
            weight[i] = 0;
            points[i].normal.x = 0;
            points[i].normal.y = 0;
        }

        // Loop over all narrow band nodes.
        for (unsigned int i=0;i<levelSet.nNarrowBand;i++)
        {
            // Node index.
            unsigned int node = levelSet.narrowBand[i];

            // Node has a neighbouring boundary point and isn't on the domain boundary.
            if ((levelSet.mesh.nodes[node].nBoundaryPoints > 0) && !levelSet.mesh.nodes[node].isDomain)
            {
                // Nodal coordinates.
                unsigned int x = levelSet.mesh.nodes[node].coord.x;
                unsigned int y = levelSet.mesh.nodes[node].coord.y;

                // The x & y gradient components.
                double gradX, gradY;

                // x direction

                // Left edge of mesh.
                if (x == 0)
                {
                    // Forward difference.
                    gradX = levelSet.signedDistance[levelSet.mesh.xyToIndex[x+1][y]]
                          - levelSet.signedDistance[levelSet.mesh.xyToIndex[x][y]];
                }

                // Right edge of mesh.
                else if (x == levelSet.mesh.width)
                {
                    // Backward difference.
                    gradX = levelSet.signedDistance[levelSet.mesh.xyToIndex[x][y]]
                          - levelSet.signedDistance[levelSet.mesh.xyToIndex[x-1][y]];
                }

                // Bulk of mesh.
                else
                {
                    // Central difference.
                    gradX = 0.5*(levelSet.signedDistance[levelSet.mesh.xyToIndex[x+1][y]]
                          - levelSet.signedDistance[levelSet.mesh.xyToIndex[x-1][y]]);
                }

                // y direction

                // Bottom edge of mesh.
                if (y == 0)
                {
                    // Forward difference.
                    gradY = levelSet.signedDistance[levelSet.mesh.xyToIndex[x][y+1]]
                          - levelSet.signedDistance[levelSet.mesh.xyToIndex[x][y]];
                }

                // Top edge of mesh.
                else if (y == levelSet.mesh.height)
                {
                    // Backward difference.
                    gradY = levelSet.signedDistance[levelSet.mesh.xyToIndex[x][y]]
                          - levelSet.signedDistance[levelSet.mesh.xyToIndex[x][y-1]];
                }

                // Bulk of mesh.
                else
                {
                    // Central difference.
                    gradY = 0.5*(levelSet.signedDistance[levelSet.mesh.xyToIndex[x][y+1]]
                          - levelSet.signedDistance[levelSet.mesh.xyToIndex[x][y-1]]);
                }

                // Absolute gradient.
                double grad = sqrt(gradX*gradX + gradY*gradY);

                // Compute the two normal vector components.
                double xNormal = gradX / grad;
                double yNormal = gradY / grad;

                // Loop over all boundary points.
                for (unsigned int j=0;j<levelSet.mesh.nodes[node].nBoundaryPoints;j++)
                {
                    // Boundary point index.
                    unsigned int point = levelSet.mesh.nodes[node].boundaryPoints[j];

                    // Distance from the boundary point to the node.
                    double dx = levelSet.mesh.nodes[node].coord.x - points[point].coord.x;
                    double dy = levelSet.mesh.nodes[node].coord.y - points[point].coord.y;

                    // Squared distance.
                    double rSqd = dx*dx + dy*dy;

                    // If boundary point lies exactly on the node, then set normal
                    // vector to that of the node.
                    if (rSqd < 1e-6)
                    {
                        points[point].normal.x = xNormal;
                        points[point].normal.y = yNormal;
                        weight[point] = 1.0;
                        isSet[point] = true;
                    }

                    else
                    {
                        // Update normal vector estimate if not already set.
                        if (!isSet[point])
                        {
                            points[point].normal.x += xNormal / rSqd;
                            points[point].normal.y += yNormal / rSqd;
                            weight[point] += 1.0 / rSqd;
                        }
                    }
                }
            }
        }

        // Compute interpolated normal vector.
        for (unsigned int i=0;i<nPoints;i++)
        {
            if (!points[i].isDomain)
            {
                points[i].normal.x /= weight[i];
                points[i].normal.y /= weight[i];

                // Compute the new vector norm.
                double norm = sqrt(points[i].normal.x*points[i].normal.x
                            + points[i].normal.y*points[i].normal.y);

                points[i].normal.x /= norm;
                points[i].normal.y /= norm;
            }
        }
    }

    double Boundary::computePerimeter(const BoundaryPoint& point)
    {
        double length = 0;

        // Sum the distance to each neighbour.
        for (unsigned int i=0;i<point.nNeighbours;i++)
        {
            double dx = point.coord.x - points[point.neighbours[i]].coord.x;
            double dy = point.coord.y - points[point.neighbours[i]].coord.y;

            length += sqrt(dx*dx + dy*dy);
        }

        return length;
    }

    void Boundary::computeMeshStatus(Mesh& mesh, const std::vector<double>* signedDistance) const
    {
        // Calculate node status.
        for (unsigned int i=0;i<mesh.nNodes;i++)
        {
            // Reset the number of boundary points associated with the element.
            mesh.nodes[i].nBoundaryPoints = 0;

            // Flag node as being on the boundary if the signed distance is within
            // a small tolerance of the zero contour. This avoids problems with
            // rounding errors when generating the discretised boundary.
            if (std::abs((*signedDistance)[i]) < 1e-6)
            {
                mesh.nodes[i].status = NodeStatus::BOUNDARY;
            }
            else if ((*signedDistance)[i] < 0)
            {
                mesh.nodes[i].status = NodeStatus::OUTSIDE;
            }
            else mesh.nodes[i].status = NodeStatus::INSIDE;
        }

        // Calculate element status.
        for (unsigned int i=0;i<mesh.nElements;i++)
        {
            // Tally counters for the element's node statistics.
            unsigned int tallyInside = 0;
            unsigned int tallyOutside = 0;

            // Reset the number of boundary segments associated with the element.
            mesh.elements[i].nBoundarySegments = 0;

            // Loop over each node of the element.
            for (unsigned int j=0;j<4;j++)
            {
                unsigned int node = mesh.elements[i].nodes[j];

                if (mesh.nodes[node].status & NodeStatus::INSIDE) tallyInside++;
                else if (mesh.nodes[node].status & NodeStatus::OUTSIDE) tallyOutside++;
            }

            // No nodes are outside: element is inside the structure.
            if (tallyOutside == 0) mesh.elements[i].status = ElementStatus::INSIDE;

            // No nodes are inside: element is outside the structure.
            else if (tallyInside == 0) mesh.elements[i].status = ElementStatus::OUTSIDE;

            // Otherwise no status.
            else mesh.elements[i].status = ElementStatus::NONE;
        }
    }

    int Boundary::isAdded(Mesh& mesh, Coord& point, const unsigned int& node,
        const unsigned int& edge, const double& distance)
    {
        // Work out the coordinates of the point (depends on which edge we are considering).
        // Set edge and distance equal to zero when considering boundary points lying exactly
        // on top of a mesh node.

        // Bottom edge.
        if (edge == 0)
        {
            point.x = mesh.nodes[node].coord.x + distance;
            point.y = mesh.nodes[node].coord.y;
        }
        // Right edge.
        else if (edge == 1)
        {
            point.x = mesh.nodes[node].coord.x;
            point.y = mesh.nodes[node].coord.y + distance;
        }
        // Top edge.
        else if (edge == 2)
        {
            point.x = mesh.nodes[node].coord.x - distance;
            point.y = mesh.nodes[node].coord.y;
        }
        // Left edge.
        else
        {
            point.x = mesh.nodes[node].coord.x;
            point.y = mesh.nodes[node].coord.y - distance;
        }

        // Check all points adjacent to the node.
        for (unsigned int i=0;i<mesh.nodes[node].nBoundaryPoints;i++)
        {
            // Index of the ith boundary point connected to the node.
            unsigned int index = mesh.nodes[node].boundaryPoints[i];

            // Point already exists.
            if ((std::abs(point.x - points[index].coord.x) < 1e-6) &&
                (std::abs(point.y - points[index].coord.y) < 1e-6))
            {
                // Boundary point is already added, return index.
                return index;
            }
        }

        // If we've made it this far, then point is new.
        return -1;
    }

    void Boundary::initialisePoint(LevelSet& levelSet, BoundaryPoint& point, const Coord& coord)
    {
        // Initialise position, length, number of segments, and isDomain.
        point.coord = coord;
        point.length = 0;
        point.nSegments = 0;
        point.nNeighbours = 0;
        point.isDomain = false;
        point.isFixed = false;

        // Assume two sensitivities to start with (objective and a single constraint).
        point.sensitivities.resize(2);

        // Initialise movement limit (CFL condition).
        point.negativeLimit = -levelSet.moveLimit;
        point.positiveLimit =  levelSet.moveLimit;

        // Check whether point lies within the move limit of the domain boundary.
        // If so, modify the lower movement limit so that point can't move outside of
        // the domain.

        // Closest distance to domain boundary in x.
        double minX = std::min(coord.x, levelSet.mesh.width - coord.x);

        // Closest distance to domain boundary in y.
        double minY = std::min(coord.y, levelSet.mesh.height - coord.y);

        // Closest distance to any domain boundary.
        double minBoundary = std::min(minX, minY);

        // Modify lower move limit.
        if (minBoundary < levelSet.moveLimit)
        {
            point.negativeLimit = -minBoundary;

            // Point is exactly on domain boundary.
            if (minBoundary < 1e-6)
                point.isDomain = true;
        }

        // Index of nearest node on the mesh.
        unsigned int node = levelSet.mesh.getClosestNode(coord);

        // Check to see if point lies on, or near, a masked region.
        if (levelSet.mesh.nodes[node].isMasked)
        {
            // Distance to masked node.
            double dx = levelSet.mesh.nodes[node].coord.x - coord.x;
            double dy = levelSet.mesh.nodes[node].coord.y - coord.y;

            // Lies on the masked node.
            if ((std::abs(dx) < 1e-6) && (std::abs(dy) < 1e-6))
                point.isFixed = true;

            // Update negative move limit.
            else
            {
                // Work out distance to the node.
                double d = sqrt(dx*dx + dy*dy);

                // Update negative move limit.
                if (-d > point.negativeLimit) point.negativeLimit = -d;
            }
        }
    }

    double Boundary::segmentLength(const BoundarySegment& segment)
    {
        // Coordinates for start and end points.
        Coord p1, p2;

        p1.x = points[segment.start].coord.x;
        p1.y = points[segment.start].coord.y;

        p2.x = points[segment.end].coord.x;
        p2.y = points[segment.end].coord.y;

        // Compute separation in x and y directions.
        double dx = p1.x - p2.x;
        double dy = p1.y - p2.y;

        return (sqrt(dx*dx + dy*dy));
    }

    void Boundary::computePointLengths()
    {
        // Loop over all boundary segments.
        for (unsigned int i=0;i<nSegments;i++)
        {
            // Add half segment length to each boundary point.
            points[segments[i].start].length += 0.5 * segments[i].length;
            points[segments[i].end].length += 0.5 * segments[i].length;

            // Update point to segment lookup for start point.
            points[segments[i].start].segments[points[segments[i].start].nSegments] = i;
            points[segments[i].start].nSegments++;

            // Update point to segment lookup for end point.
            points[segments[i].end].segments[points[segments[i].end].nSegments] = i;
            points[segments[i].end].nSegments++;

            // Update nearest neighbours for the start point.
            points[segments[i].start].neighbours[points[segments[i].start].nNeighbours] = segments[i].end;
            points[segments[i].start].nNeighbours++;

            // Update nearest neighbours for the end point.
            points[segments[i].end].neighbours[points[segments[i].end].nNeighbours] = segments[i].start;
            points[segments[i].end].nNeighbours++;
        }
    }
}
