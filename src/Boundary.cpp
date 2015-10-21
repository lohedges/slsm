/*
 * File:	Boundary.cpp
 * Author:	lester
 */

#include "Boundary.h"

Boundary::Boundary(Mesh& mesh_, LevelSet& levelSet_) : mesh(mesh_), levelSet(levelSet_)
{
}

void Boundary::discretise()
{
    // Reserve memory for boundary points and segments.
    points.reserve(0.2*mesh.nNodes);
    segments.reserve(0.2*mesh.nNodes);

    // Reset the number of points and segments.
    nPoints = nSegments = 0;

    // Compute the status of nodes and elements in the finite element mesh.
    computeMeshStatus();

    // Loop over all elements.
    for (unsigned int i=0;i<mesh.nElements;i++)
    {
        // The element isn't outside of the structure.
        if (mesh.elements[i].status != ElementStatus::OUTSIDE)
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
                unsigned int n1 = mesh.elements[i].nodes[j];

                // Index of second node (reconnecting to 0th node).
                // Edge connectivity goes: 0 --> 1, 1 --> 2, 2 --> 3, 3 --> 0
                unsigned int n2 = (j == 3) ? 0 : (j + 1);

                // Conver to node index.
                n2 = mesh.elements[i].nodes[n2];

                // One node is inside, the other is outside. The edge is cut.
                if ((mesh.nodes[n1].status|mesh.nodes[n2].status) == NodeStatus::CUT)
                {
                    // Compute the distance from node 1 to the intersection point (by interpolation).
                    double d = levelSet.signedDistance[n1]
                             / (levelSet.signedDistance[n1] - levelSet.signedDistance[n2]);

                    d = std::abs(d);

                    // Create the boundary point.
                    Coord point;

                    // Work out the coordinates of the point (depends on which edge we are considering).
                    if (j == 0)
                    {
                        point.x = mesh.nodes[n1].coord.x + d;
                        point.y = mesh.nodes[n1].coord.y;
                    }
                    else if (j == 1)
                    {
                        point.x = mesh.nodes[n1].coord.x;
                        point.y = mesh.nodes[n1].coord.y + d;
                    }
                    else if (j == 2)
                    {
                        point.x = mesh.nodes[n1].coord.x - d;
                        point.y = mesh.nodes[n1].coord.y;
                    }
                    else
                    {
                        point.x = mesh.nodes[n1].coord.x;
                        point.y = mesh.nodes[n1].coord.y - d;
                    }

                    // Make sure that the boundary point hasn't already been added.
                    bool isAdded = false;

                    for (unsigned int k=0;k<mesh.nodes[n1].nBoundaryPoints;k++)
                    {
                        // Index of the kth boundary point connected to node n1.
                        unsigned int index = mesh.nodes[n1].boundaryPoints[k];

                        // Point already exists.
                        if ((std::abs(point.x - points[index].x) < 1e-6) &&
                            (std::abs(point.y - points[index].y) < 1e-6))
                        {
                            // Store existing boundary point.
                            boundaryPoints[nCut] = index;

                            // Flag boundary point as already added.
                            isAdded = true;
                            break;
                        }
                    }

                    if (!isAdded)
                    {
                        // Add boundary point.
                        points.push_back(point);

                        // Create node to boundary point lookup.
                        mesh.nodes[n1].boundaryPoints[mesh.nodes[n1].nBoundaryPoints] = nPoints;
                        mesh.nodes[n2].boundaryPoints[mesh.nodes[n2].nBoundaryPoints] = nPoints;
                        mesh.nodes[n1].nBoundaryPoints++;
                        mesh.nodes[n2].nBoundaryPoints++;

                        // Store boundary point for cut edge.
                        boundaryPoints[nCut] = nPoints;

                        // Increment number of boundary points.
                        nPoints++;
                    }

                    // Increment number of cut edges.
                    nCut++;
                }

                // Both nodes lie on the boundary.
                else if ((mesh.nodes[n1].status & NodeStatus::BOUNDARY) &&
                    (mesh.nodes[n2].status & NodeStatus::BOUNDARY))
                {
                    // Create boundary segment.
                    BoundarySegment segment;
                    segment.node1 = n1;
                    segment.node2 = n2;
                    segment.element = i;

                    // Create element to segment lookup.
                    mesh.elements[i].boundarySegments[mesh.elements[i].nBoundarySegments] = nSegments;
                    mesh.elements[i].nBoundarySegments++;

                    // Add segment to vector.
                    segments.push_back(segment);
                    nSegments++;
                }
            }

            // If the element was cut, determine the boundary segment(s).

            // If two edges are cut, then a boundary segment must cross both.
            if (nCut == 2)
            {
                // Create boundary segment.
                // Here the start and end "nodes" are actually boundary points.
                BoundarySegment segment;
                segment.node1 = mesh.nNodes + boundaryPoints[0];
                segment.node2 = mesh.nNodes + boundaryPoints[1];
                segment.element = i;

                // Create element to segment lookup.
                mesh.elements[i].boundarySegments[mesh.elements[i].nBoundarySegments] = nSegments;
                mesh.elements[i].nBoundarySegments++;

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
                    unsigned int node = mesh.elements[i].nodes[j];

                    // Node is on the boundary, check its neighbours.
                    if (mesh.nodes[node].status & NodeStatus::BOUNDARY)
                    {
                        // Index of next node.
                        unsigned int nAfter = (j == 3) ? 0 : (j + 1);

                        // Index of previous node.
                        unsigned int nBefore = (j == 0) ? 3 : (j - 1);

                        // Convert to node indices.
                        nAfter = mesh.elements[i].nodes[nAfter];
                        nBefore = mesh.elements[i].nodes[nBefore];

                        // If a neighbour is outside the boundary, then add a boundary segment.
                        if ((mesh.nodes[nAfter].status & NodeStatus::OUTSIDE) ||
                            (mesh.nodes[nBefore].status & NodeStatus::OUTSIDE))
                        {
                            // Create boundary segment.
                            BoundarySegment segment;
                            segment.node1 = mesh.nNodes + boundaryPoints[0];
                            segment.node2 = node;
                            segment.element = i;

                            // Create element to segment lookup.
                            mesh.elements[i].boundarySegments[mesh.elements[i].nBoundarySegments] = nSegments;
                            mesh.elements[i].nBoundarySegments++;

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
                    unsigned int node = mesh.elements[i].nodes[j];

                    lsfSum += levelSet.signedDistance[node];
                }

                // Create boundary segment.
                BoundarySegment segment;

                // Store the status of the first node.
                unsigned int node = mesh.elements[i].nodes[0];
                NodeStatus::NodeStatus status = mesh.nodes[node].status;

                if (((status & NodeStatus::INSIDE) && (lsfSum > 0)) ||
                    ((status & NodeStatus::OUTSIDE) && (lsfSum < 0)))
                {
                    segment.node1 = mesh.nNodes + boundaryPoints[0];
                    segment.node2 = mesh.nNodes + boundaryPoints[1];
                    segment.element = i;

                    // Create element to segment lookup.
                    mesh.elements[i].boundarySegments[mesh.elements[i].nBoundarySegments] = nSegments;
                    mesh.elements[i].nBoundarySegments++;

                    // Add segment to vector.
                    segments.push_back(segment);
                    nSegments++;

                    segment.node1 = mesh.nNodes + boundaryPoints[2];
                    segment.node2 = mesh.nNodes + boundaryPoints[3];
                    segment.element = i;

                    // Create element to segment lookup.
                    mesh.elements[i].boundarySegments[mesh.elements[i].nBoundarySegments] = nSegments;
                    mesh.elements[i].nBoundarySegments++;

                    // Add segment to vector.
                    segments.push_back(segment);
                    nSegments++;
                }

                else
                {
                    segment.node1 = mesh.nNodes + boundaryPoints[0];
                    segment.node2 = mesh.nNodes + boundaryPoints[3];
                    segment.element = i;

                    // Create element to segment lookup.
                    mesh.elements[i].boundarySegments[mesh.elements[i].nBoundarySegments] = nSegments;
                    mesh.elements[i].nBoundarySegments++;

                    // Add segment to vector.
                    segments.push_back(segment);
                    nSegments++;

                    segment.node1 = mesh.nNodes + boundaryPoints[1];
                    segment.node2 = mesh.nNodes + boundaryPoints[2];
                    segment.element = i;

                    // Create element to segment lookup.
                    mesh.elements[i].boundarySegments[mesh.elements[i].nBoundarySegments] = nSegments;
                    mesh.elements[i].nBoundarySegments++;

                    // Add segment to vector.
                    segments.push_back(segment);
                    nSegments++;
                }

                // Update element status to indicate whether centre is in or out.
                mesh.elements[i].status = (lsfSum > 0) ? ElementStatus::CENTRE_INSIDE : ElementStatus::CENTRE_OUTSIDE;
            }

            // If no edges are cut and element is not inside structure
            // then the boundary segment must cross the diagonal.
            else if ((nCut == 0) && (mesh.elements[i].status != ElementStatus::INSIDE))
            {
                // Find the two boundary nodes.
                for (unsigned int j=0;j<4;j++)
                {
                    // Node index.
                    unsigned int node = mesh.elements[i].nodes[j];

                    if (mesh.nodes[node].status & NodeStatus::BOUNDARY)
                    {
                        boundaryPoints[nCut] = node;
                        nCut++;
                    }
                }

                // Create boundary segment.
                BoundarySegment segment;
                segment.node1 = boundaryPoints[0];
                segment.node2 = boundaryPoints[1];
                segment.element = i;

                // Create element to segment lookup.
                mesh.elements[i].boundarySegments[mesh.elements[i].nBoundarySegments] = nSegments;
                mesh.elements[i].nBoundarySegments++;

                // Add segment to vector.
                segments.push_back(segment);
                nSegments++;
            }
        }
    }
}

void Boundary::computeMeshStatus()
{
    // Calculate node status.
    for (unsigned int i=0;i<mesh.nNodes;i++)
    {
        // Reset the number of boundary points associated with the element.
        mesh.nodes[i].nBoundaryPoints = 0;

        if (levelSet.signedDistance[i] == 0)
        {
            mesh.nodes[i].status = NodeStatus::BOUNDARY;
        }
        else if (levelSet.signedDistance[i] < 0)
        {
            mesh.nodes[i].status = NodeStatus::INSIDE;
        }
        else mesh.nodes[i].status = NodeStatus::OUTSIDE;
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
