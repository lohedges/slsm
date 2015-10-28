/*
 * File:    Boundary.h
 * Author:  lester
 */

#ifndef _BOUNDARY_H
#define _BOUNDARY_H

#include <functional>

#include "Common.h"
#include "LevelSet.h"

/*! \file Boundary.h
    \brief A class computing and storing the discretised boundary.
 */

// ASSOCIATED DATA TYPES

//! \brief A container for storing information associated with a boundary segment.
struct BoundarySegment
{
    unsigned int start;         //!< Index of start point (node or boundary point).
    unsigned int end;           //!< Index of end point (node or boundary point).
    unsigned int element;       //!< The element cut by the boundary segment.
    double length;              //!< Length of the boundary segment.
    double weight;              //!< Weighting factor for boundary segment.
};

// MAIN CLASS

/*! A class computing and storing the discretised boundary.

    The boundary is computed by looking for nodes lying exactly on the zero
    contour of the level set and then constructing a set of additional boundary
    points by simple linear interpolation when the level set changes sign
    between the nodes on an element edge.

    The points vector holds coordinates for the interpolated boundary points.
    Boundary segment data is stored in the segments vector. Note that the
    start and end points of a boundary segment can either be a node, or a
    boundary point. If a point is a node, the index will lie in the range
    0 to mesh.nNodes - 1, if it's a boundary point then it will be between
    mesh.nNodes and mesh.nNodes + nPoints - 1.
 */
class Boundary
{
public:
    //! Constructor.
    /*! \param mesh_
            A reference to the finite element mesh.

        \param levelSet_
            A reference to the level set object.
     */
    Boundary(Mesh&, LevelSet&);

    //! Use linear interpolation to compute the discretised boundary
    //! (zero contour) of the level set.
    void discretise();

    //! Calculate the material area fraction in each element.
    void computeAreaFractions();

    /// Vector of boundary points.
    std::vector<Coord> points;

    /// Vector of boundary segments.
    std::vector<BoundarySegment> segments;

    /// The number of boundary points.
    unsigned int nPoints;

    /// The number of boundary segments.
    unsigned int nSegments;

    /// The total length of the boundary.
    double length;

    /// The total area fraction of the mesh.
    double area;

private:
    /// A reference to the finite element mesh.
    Mesh& mesh;

    /// A reference to the level set object.
    LevelSet& levelSet;

    //! Determine the status of the elements and nodes of the
    //! finite element grid.
    void computeMeshStatus();

    //! Calculate the material area for an element cut by the boundary.
    /*! \param element
            A reference to the element.

        \return
            The area fraction.
     */
    double cutArea(const Element&);

    //! Whether a point is clockwise of another.
    /*! \param point1
            The coordinates of the first point.

        \param point2
            The coordinates of the second point.

        \param centre
            The coordinates of the element centre.

        \return
            Whether the first point is clockwise of the second.
     */
    bool isClockwise(const Coord&, const Coord&, const Coord&) const;

    //! Return the area of a polygon.
    /*! \param vertices
            A clockwise ordered vector of polygon vertices.

        \param nVertices
            The number of vertices.

        \param centre
            The coordinates of the element centre.

        \return
            The area of the polygon.
     */
    double polygonArea(std::vector<Coord>&, const unsigned int&, const Coord&) const;

    //! Return the length of a boundary segment.
    /*! \param segment
            A reference to the boundary segment.

        \return
            The length of the bounary segment.
     */
    double segmentLength(const BoundarySegment&);
};

#endif  /* _BOUNDARY_H */
