/*
 * File:	Boundary.h
 * Author:	lester
 */

#ifndef _BOUNDARY_H
#define _BOUNDARY_H

#include "Common.h"
#include "LevelSet.h"

struct BoundarySegment
{
    unsigned int node1;         //!< Index for first node end point.
    unsigned int node2;         //!< Index for second node end point.
    unsigned int element;       //!< The element cut by the boundary segment.
    double length;              //!< Length of the boundary segment.
    double weight;              //!< Weighting factor for boundary segment.
};

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

    /// Vector of boundary points.
    std::vector<Coord> points;

    /// Vector of boundary segments.
    std::vector<BoundarySegment> segments;

    /// The number of boundary points.
    unsigned int nPoints;

    /// The number of boundary segments.
    unsigned int nSegments;

private:
    /// A reference to the finite element mesh.
    Mesh& mesh;

    /// A reference to the level set object.
    LevelSet& levelSet;

    //! Determine the status of the elements and nodes of the
    //! finite element grid.
    void computeMeshStatus();
};

#endif	/* _BOUNDARY_H */
