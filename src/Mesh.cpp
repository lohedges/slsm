/*
 * File:    Mesh.cpp
 * Author:  lester
 */

#include "Mesh.h"

Mesh::Mesh(unsigned int width_,
           unsigned int height_,
           bool isPeriodic_) :

           width(width_),
           height(height_),
           nElements(width*height),
           nNodes((1+width)*(1+height)),
           isPeriodic(isPeriodic_)
{
    // Resize element and node data structures.
    elements.resize(nElements);
    nodes.resize(nNodes);

    // Calculate node nearest neighbours.
    initialiseNodes();

    // Initialise elements (and node to element connectivity).
    initialiseElements();
}

unsigned int Mesh::getClosestNode(double x, double y)
{
    return 0;
}

unsigned int Mesh::getElement(double x, double y)
{
    return 0;
}

void Mesh::initialiseNodes()
{
    // Coordinates of the node.
    unsigned int x, y;

    // Loop over all nodes.
    for (unsigned int i=0;i<nNodes;i++)
    {
        // Zero number of connected elements.
        nodes[i].nElements = 0;

        // Work out node coordinates.
        x = i % (width + 1);
        y = int(i / (width + 1));

        // Set node coordinates.
        nodes[i].coord.x = x;
        nodes[i].coord.y = y;

        // Determine nearest neighbours.
        initialiseNeighbours(i, x, y);
    }
}

void Mesh::initialiseElements()
{
    // Coordinates of the element.
    unsigned int x, y;

    // Number of nodes along width of mesh (number of elements plus one)
    unsigned int w = width + 1;

    // Loop over all elements.
    for (unsigned int i=0;i<nElements;i++)
    {
        // Work out element coordinates.
        x = i % width;
        y = int(i / width);

        // Store coordinates of elemente centre.
        elements[i].coord.x = x + 0.5;
        elements[i].coord.y = y + 0.5;

        // Store connectivity (element --> node)

        // Node on bottom left corner of element.
        elements[i].nodes[0] = x + (y * w);

        // Node on bottom right corner of element.
        elements[i].nodes[1] = x + 1 + (y * w);

        // Node on top right corner of element.
        elements[i].nodes[2] = x + 1 + ((y + 1) * w);

        // Node on top right corner of element.
        elements[i].nodes[3] = x + ((y + 1) * w);

        // Fill reverse connectivity arrays (node --> element)
        for (unsigned int j=0;j<4;j++)
        {
            unsigned int node = elements[i].nodes[j];
            nodes[node].elements[nodes[node].nElements] = i;
            nodes[node].nElements++;
        }
    }
}

void Mesh::initialiseNeighbours(unsigned int node, unsigned int x, unsigned int y)
{
    // Number of nodes along width of mesh (number of elements plus one)
    unsigned int w = width + 1;

    // Neighbours to left and right.
    nodes[node].neighbours[0] = (x - 1 + w) % w + (y * w);
    nodes[node].neighbours[1] = (x + 1 + w) % w + (y * w);

    // Neighbours below and above.
    nodes[node].neighbours[2] = x + (w * ((y - 1 + w) % w));
    nodes[node].neighbours[3] = x + (w * ((y + 1 + w) % w));

    // The mesh isn't periodic, flag out of bounds neigbours.
    if (!isPeriodic)
    {
        // Node is on first or last row.
        if (x == 0) nodes[node].neighbours[0] = nNodes;
        else if (x == width) nodes[node].neighbours[1] = nNodes;

        // Node is on first or last column.
        if (y == 0) nodes[node].neighbours[2] = nNodes;
        else if (y == height) nodes[node].neighbours[3] = nNodes;
    }
}
