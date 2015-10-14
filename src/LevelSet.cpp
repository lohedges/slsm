/*
 * File:	LevelSet.cpp
 * Author:	lester
 */

#include "LevelSet.h"

LevelSet::LevelSet(Mesh& mesh, unsigned int bandWidth_) :
    nNodes(mesh.nNodes),
    bandWidth(bandWidth_)
{
    errno = EINVAL;
    check(bandWidth > 2, "Width of the narrow band must be greater than 2.");

    // Resize data structures.
    signedDistance.resize(nNodes);
    gradient.resize(nNodes);
    narrowBand.resize(nNodes);

    // Resize mines vector, 20% of nNodes is a reasonable estimate.
    mines.resize(int(0.2*nNodes));

    // Generate a Swiss cheese structure.
    initialise(mesh);

    // Initialise the narrow band.
    initialiseNarrowBand(mesh);

    return;

error:
    exit(EXIT_FAILURE);
}

LevelSet::LevelSet(Mesh& mesh, unsigned int bandWidth_, const std::vector<Hole>& holes) :
    nNodes(mesh.nNodes),
    bandWidth(bandWidth_)
{
    errno = EINVAL;
    check(bandWidth > 2, "Width of the narrow band must be greater than 2.");

    // Resize data structures.
    signedDistance.resize(nNodes);
    gradient.resize(nNodes);
    narrowBand.resize(nNodes);

    // Resize mines vector, 20% of nNodes is a reasonable estimate.
    mines.resize(int(0.2*nNodes));

    // Initialise level set function from hole array.
    initialise(mesh, holes);

    // Initialise the narrow band.
    initialiseNarrowBand(mesh);

    return;

error:
    exit(EXIT_FAILURE);
}

void LevelSet::update()
{
}

void LevelSet::initialise(const Mesh& mesh)
{
    // First initialise LSF based on domain boundary.
    closestDomainBoundary(mesh);

    // Now add Swiss cheese hole arrangement.
}

void LevelSet::initialise(const Mesh& mesh, const std::vector<Hole>& holes)
{
    // First initialise LSF based on domain boundary.
    closestDomainBoundary(mesh);

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

void LevelSet::closestDomainBoundary(const Mesh& mesh)
{
    // Initial LSF is distance from closest domain boundary.
    for (unsigned int i=0;i<nNodes;i++)
    {
        // Closest edge in x.
        unsigned int minX = std::min(mesh.nodes[i].coord.x, mesh.width - mesh.nodes[i].coord.x);

        // Closest edge in y.
        unsigned int minY = std::min(mesh.nodes[i].coord.y, mesh.width - mesh.nodes[i].coord.y);

        // Signed distance is the minimum of minX and minY;
        signedDistance[i] = double(std::min(minX, minY));
    }
}

void LevelSet::initialiseNarrowBand(Mesh& mesh)
{
    unsigned int mineWidth = bandWidth - 1;

    // Reset the number of nodes in the narrow band.
    nNarrowBand = 0;

    // Reset the number of mines.
    nMines = 0;

    // Reset minimum distance to zero isocontour.
    minDistance = 1e6;

    // Loop over all nodes.
    for (unsigned int i=0;i<nNodes;i++)
    {
        // Absolute value of the signed distance function.
        double absoluteSignedDistance = std::abs(signedDistance[i]);

        // Check whether the current node is the closest
        // to the boundary (zero isocontour).
        if (absoluteSignedDistance < minDistance)
        {
            minDistance = absoluteSignedDistance;
            boundaryNode = i;
        }

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
