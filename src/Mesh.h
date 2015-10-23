/*
 * File:	Mesh.h
 * Author:	lester
 */

#ifndef _MESH_H
#define _MESH_H

#include <vector>

#include "Common.h"

/*! \file Mesh.h
    \brief A class for the fixed grid, finite-element mesh.
 */

// ASSOCIATED DATA TYPES

//! Whether a node lies inside, outside, or on the boundary.
namespace NodeStatus
{
    // Left bit shift enumerated types to allow the creation
    // of sets and simple bit masking operations.
    // For example, to test whether a node is neither inside or outside...
    // if (node.status ! (NodeStatus::Inside|NodeStatus::Outside))
    enum NodeStatus
    {
        NONE            = 0,                //!< No status.
        INSIDE          = (1 << 0),         //!< Node lies inside the boundary.
        OUTSIDE         = (1 << 1),         //!< Node lies outside the boundary.
        BOUNDARY        = (1 << 2),         //!< Node is cut by the boundary.
        MASKED          = (1 << 3),         //!< Node is masked (for inserting obstacles, etc.)
        CUT             = (INSIDE|OUTSIDE), //!< Node pair cuts the boundary.
    };
}

//! Whether a element lies completely inside or outside the structure.
//! If not, whether the element centre lies inside or outside.
namespace ElementStatus
{
    // Left bit shift enumerated types to allow the creation
    // of sets and simple bit masking operations.
    enum ElementStatus
    {
        NONE            = 0,                //!< No status.
        INSIDE          = (1 << 0),         //!< Element is entirely inside the structure.
        OUTSIDE         = (1 << 1),         //!< Element is entirely outside the structure.
        CENTRE_INSIDE   = (1 << 2),         //!< Element centre lies inside the structure.
        CENTRE_OUTSIDE  = (1 << 3),         //!< Element centre lies outside the structure.
    };
}

//! Structure containing attributes for an individual grid element.
struct Element
{
    double area;                            //!< Material area fraction.
    unsigned int nodes[4];                  //!< Indices for nodes of the element.
    unsigned int boundarySegments[2];       //!< Indices for boundary segments associated with the element.
    unsigned int nBoundarySegments;         //!< The number of boundary segments associated with the element.
    ElementStatus::ElementStatus status;    //!< Whether the element (or it's centre) lies inside or outside the structure.
};

//! Structure containing attributes for an individual grid node.
struct Node
{
    Coord coord;                            //!< Node coordinate.
    unsigned int neighbours[4];             //!< Indices of nearest neighbour nodes.
    unsigned int elements[4];               //!< Indices of elements the node is connected to.
    unsigned int nElements;                 //!< Number of elements that the node is connected to (no periodicity).
    unsigned int boundaryPoints[4];         //!< Indices of boundary points associated with the node.
    unsigned int nBoundaryPoints;           //!< The number of boundary points associated with the node.
    bool isActive;                          //!< Whether the node is active (part of narrow band, and not fixed).
    bool isMine;                            //!< Whether the node lies on the edge of the narrow band.
    NodeStatus::NodeStatus status;          //!< Whether node is outside, inside, or on the boundary.
};

// MAIN CLASS

/*!\brief A class for the fixed grid, finite-element mesh.

    Stores connectivty information between grid elements and nodes.
    The grid is assumed to be two-dimensional and is comprised of
    square elements of unit side.

    Elements are comprised of four nodes, labelled in anticlockwise order
    from the bottom left, i.e. bottom left, bottom right, top right, top left.

    Each node has four nearest neighbours ordered as left, right, down, up.
    Diagonal neighbours can be accessed by looking at neighbours of neighbours,
    e.g. for the lower left diagonal of node i

    node = nodes[nodes[i].neighbours[0]].neigbours[2];

    Note that diagonal nodes can be accessed in multiple ways, e.g.
    for the lower left diagonal of node i we could go left then down,
    or down then left (shortest paths).

    For non-periodic meshes, neighbours that are outside of the domain are
    given the value nNodes, i.e. one past the end of the node array, which
    runs from 0 to nNodes - 1.
 */
class Mesh
{
public:
    //! Constructor.
    /*! \param width_
            The width of the finite element grid.

        \param height_
            The height of the finite element grid.

        \param isPeriodic_
            Whether the grid is periodic.
     */
	Mesh(unsigned int, unsigned int, bool);

    //! For a given x-y coordinate, find the index of the closest node.
    /*! \param x
            The x position of the point.

        \param y
            The y position of the point.

        \return
            The index of the closest node.
     */
    unsigned int getClosestNode(double, double);

    //! For a given x-y coordinate, find the element that contains the point.
    /*! \param x
            The x position of the point.

        \param y
            The y position of the point.

        \return
            The index of containing element.
     */
    unsigned int getElement(double, double);

    std::vector<Element> elements;  //!< Fixed grid elements (cells).
    std::vector<Node> nodes;        //!< Fixed grid nodes.

    const unsigned int width;       //!< The grid width (number of elements in x).
    const unsigned int height;      //!< The grid height (number of elements in y).
    const unsigned int nElements;   //!< The total number of grid elements.
    const unsigned int nNodes;      //!< The total number of nodes.

private:
    const bool isPeriodic;          //!< Whether the grid has periodic boundaries.
    const double tolerance = 1e-6;  //!< Rounding tolerance.

    //! Initialise mesh nodes.
    void initialiseNodes();

    //! Initialise mesh elements.
    void initialiseElements();

    //! Compute the nearest neighbours of a node.
    /*! \param node
            The node index.

        \param x
            The x coordinate of the node.

        \param y
            The y coordinate of the node.
     */
    void initialiseNeighbours(unsigned int, unsigned int, unsigned int);
};

#endif	/* _MESH_H */
