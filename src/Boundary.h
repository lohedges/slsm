/*
  Copyright (c) 2015-2016 Lester Hedges <lester.hedges+lsm@gmail.com>

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

#ifndef _BOUNDARY_H
#define _BOUNDARY_H

#include <algorithm>
#include <functional>

#include "Common.h"
#include "LevelSet.h"

/*! \file Boundary.h
    \brief A class for the discretised boundary.
 */

namespace lsm
{
    // ASSOCIATED DATA TYPES

    //! \brief A container for storing information associated with a boundary point.
    struct BoundaryPoint
    {
        Coord coord;                        //!< Coordinate of the boundary point.
        double length;                      //!< Integral length of the boundary point.
        double curvature;                   //!< Local interface curvature around point.
        double velocity;                    //!< Normal velocity (positive acts inwards).
        double negativeLimit;               //!< Movement limit in negative direction (inwards).
        double positiveLimit;               //!< Movement limit in positive direction (outwards).
        bool isDomain;                      //!< Whether the point lies close to the domain boundary.
        unsigned int nSegments;             //!< The number of boundary segments that a point belongs to.
        unsigned int segments[2];           //!< The indices of the two segments to which a point belongs.
        std::vector<double> sensitivities;  //!< Objective and constraint sensitivities.
    };

    //! \brief A container for storing information associated with a boundary segment.
    struct BoundarySegment
    {
        unsigned int start;                 //!< Index of start point.
        unsigned int end;                   //!< Index of end point.
        unsigned int element;               //!< The element cut by the boundary segment.
        double length;                      //!< Length of the boundary segment.
        double weight;                      //!< Weighting factor for boundary segment.
    };

    // MAIN CLASS

    /*! \brief A class computing and storing the discretised boundary.

        The boundary is computed by looking for nodes lying exactly on the zero
        contour of the level set and then constructing a set of additional boundary
        points by simple linear interpolation when the level set changes sign
        between the nodes on an element edge.

        The points vector holds coordinates for boundary points (both those lying
        exactly on nodes of the finite element mesh, and the interpolated points).
        Boundary segment data is stored in the segments vector.
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
        /*! \param isTarget
                Whether to discretise the target signed distance function.
         */
        void discretise(bool isTarget = false);

        //! Calculate the material area fraction in each element.
        /*! \return
                The total element area fraction.
         */
        double computeAreaFractions();

        //! Compute the local curvature around each boundary point.
        /*! \return
                The mean boundary point curvature.
         */
        double computeCurvatures();

        /// Vector of boundary points.
        std::vector<BoundaryPoint> points;

        /// Vector of boundary segments.
        std::vector<BoundarySegment> segments;

        /// The number of boundary points.
        unsigned int nPoints;

        /// The number of boundary segments.
        unsigned int nSegments;

        /// The total length of the boundary.
        double length;

        /// The mean curvature.
        double curvature;

        /// The total area fraction of the mesh.
        double area;

    private:
        /// A reference to the finite element mesh.
        Mesh& mesh;

        /// A reference to the level set object.
        LevelSet& levelSet;

        //! Determine the status of the elements and nodes of the finite element grid.
        /*! \param signedDistance
                A pointer to the signed distance function vector.
         */
        void computeMeshStatus(const std::vector<double>* signedDistance) const;

        //! Check whether a boundary point has already been added.
        /*! \param point
                The coordinates of the boundary point (to be determined).

            \param node
                The index of the adjacent node.

            \param edge
                The index of the element edge.

            \param distance
                The distance from the node.

            \return
                The index of the boundary point if previously added, minus one if not.
         */
        int isAdded(Coord&, const unsigned int&, const unsigned int&, const double&);

        //! Initialise a boundary point.
        /*! \param point
                A reference to a boundary point.

            \param coord
                The position vector of the boundary point.
         */
        void initialisePoint(BoundaryPoint&, const Coord&);

        //! Calculate the material area for an element cut by the boundary.
        /*! \param element
                A reference to the element.

            \return
                The area fraction.
         */
        double cutArea(const Element&);

        //! Whether a point is clockwise of another. The origin point is 12 o'clock.
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

        //! Compute the (potentially weighted) integral length for each boundary point.
        void computePointLengths();
    };
}

#endif  /* _BOUNDARY_H */
