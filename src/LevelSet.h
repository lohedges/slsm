/*
  Copyright (c) 2015-2017 Lester Hedges <lester.hedges+slsm@gmail.com>

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

#ifndef _LEVELSET_H
#define _LEVELSET_H

#include "Common.h"
#include "Mesh.h"

/*! \file LevelSet.h
    \brief A class for the level set function.
 */

namespace slsm
{
    // FORWARD DECLARATIONS

    class  Boundary;
    struct BoundaryPoint;
    class  Hole;
    class  MersenneTwister;

    /*! \brief A class for the level set function.

        The level set is represented as a signed distance function from the zero
        countour. Positive values are inside the structure, negative values are
        outside.

        The class provides methods for initialising and reinitialising the signed
        distance function. The default initialisation uses a set of linearly spaced
        circular holes to create the classic "Swiss cheese" starting configuration.
        Alternatively, the user can pass the constructor a vector of holes or point
        coordinates that will be used to define the initial structure. A vector
        of coordinates defining the interface of a target structure can be used
        for shape matching simulations.

        Reinitialisation of the signed distance function is performed using an
        implementation of the Fast Marching Method. A second-order stencil is used
        for the upwind finite difference scheme where possible.

        Functionality is also provided for tracking nodes that are part of the
        narrow band region around the zero contour, as well as mine nodes at
        the edge of the narrow band.
     */
    class LevelSet
    {
    public:
        //! Constructor.
        /*! \param width
                The width of the fixed-grid mesh.

            \param height
                The height of the fixed-grid mesh.

            \param moveLimit_
                The CFL limit (in units of the mesh grid spacing).

            \param bandWidth_
                The width of the narrow band region.

            \param isFixedDomain_
                Whether the domain boundary is fixed.
         */
        LevelSet(unsigned int, unsigned int, double moveLimit_ = 0.5,
            unsigned int bandWidth_ = 6, bool isFixedDomain_ = false);

        //! Constructor.
        /*! \param width
                The width of the fixed-grid mesh.

            \param height
                The height of the fixed-grid mesh.

            \param holes
                A vector of holes.

            \param moveLimit_
                The CFL limit (in units of the mesh grid spacing).

            \param bandWidth_
                The width of the narrow band region.

            \param isFixedDomain_
                Whether the domain boundary is fixed.
         */
        LevelSet(unsigned int, unsigned int, const std::vector<Hole>&, double moveLimit_ = 0.5,
            unsigned int bandWidth_ = 6, bool isFixedDomain_ = false);

        //! Constructor.
        /*! \param width
                The width of the fixed-grid mesh.

            \param height
                The height of the fixed-grid mesh.

            \param points
                A vector of point coordinates (clockwise ordered and closed).

            \param moveLimit_
                The CFL limit (in units of the mesh grid spacing).

            \param bandWidth_
                The width of the narrow band region.

            \param isFixedDomain_
                Whether the domain boundary is fixed.
         */
        LevelSet(unsigned int, unsigned int, const std::vector<Coord>&, double moveLimit_ = 0.5,
            unsigned int bandWidth_ = 6, bool isFixedDomain_ = false);

        //! Constructor.
        /*! \param width
                The width of the fixed-grid mesh.

            \param height
                The height of the fixed-grid mesh.

            \param initialHoles
                A vector of holes for the initial interface.

            \param targetHoles
                A vector of holes for the target interface.

            \param moveLimit_
                The CFL limit (in units of the mesh grid spacing).

            \param bandWidth_
                The width of the narrow band region.

            \param isFixedDomain_
                Whether the domain boundary is fixed.
         */
        LevelSet(unsigned int, unsigned int, const std::vector<Hole>&, const std::vector<Hole>&,
            double moveLimit_ = 0.5, unsigned int bandWidth_ = 6, bool isFixedDomain_ = false);

        //! Constructor.
        /*! \param width
                The width of the fixed-grid mesh.

            \param height
                The height of the fixed-grid mesh.

            \param holes
                A vector of holes for the initial interface.

            \param points
                A vector of point coordinates for the target interface (clockwise ordered and closed).

            \param moveLimit_
                The CFL limit (in units of the mesh grid spacing).

            \param bandWidth_
                The width of the narrow band region.

            \param isFixedDomain_
                Whether the domain boundary is fixed.
         */
        LevelSet(unsigned int, unsigned int, const std::vector<Hole>&, const std::vector<Coord>&,
            double moveLimit_ = 0.5, unsigned int bandWidth_ = 6, bool isFixedDomain_ = false);

        //! Constructor.
        /*! \param width
                The width of the fixed-grid mesh.

            \param height
                The height of the fixed-grid mesh.

            \param initialPoints
                A vector of point coordinates for the initial interface (clockwise ordered and closed).

            \param targetPoints
                A vector of point coordinates for the target interface (clockwise ordered and closed).

            \param moveLimit_
                The CFL limit (in units of the mesh grid spacing).

            \param bandWidth_
                The width of the narrow band region.

            \param isFixedDomain_
                Whether the domain boundary is fixed.
         */
        LevelSet(unsigned int, unsigned int, const std::vector<Coord>&, const std::vector<Coord>&,
            double moveLimit_ = 0.5, unsigned int bandWidth_ = 6, bool isFixedDomain_ = false);

        //! Update the level set function.
        /*! \param timeStep
                The time step.

            \return
                Whether the signed distance was reinitialised.
         */
        bool update(double);

        //! Mask off a region of the domain.
        /*! param holes
                A reference to a vector of holes.
         */
        void mask(const std::vector<Hole>&);

        //! Mask off a region of the domain.
        /*! param points
                A reference to a vector of points.
         */
        void mask(const std::vector<Coord>&);

        //! Reinitialise the level set to a signed distance function.
        void reinitialise();

        //! Extend boundary point velocities to the level set nodes.
        /*! \param boundaryPoints
                A reference to a vector of boundary points.
         */
        void computeVelocities(const std::vector<BoundaryPoint>&);

        //! Extend boundary point velocities to the level set nodes.
        /*! \param boundaryPoints
                A reference to a vector of boundary points.

            \param timeStep
                The time step for the level set update.

            \param temperature
                The temperature of the thermal bath.

            \param rng
                A reference to the random number generator.

            \return
                The time step scaling factor.
         */
        double computeVelocities(std::vector<BoundaryPoint>&, double&, const double, MersenneTwister&);

#ifdef PYBIND
        //! Extend boundary point velocities to the level set nodes.
        /*! \param boundaryPoints
                A reference to a vector of boundary points.

            \param timeStep
                The time step for the level set update.

            \param temperature
                The temperature of the thermal bath.

            \param rng
                A reference to the random number generator.

            \return
                The time step scaling factor.
         */
        double computeVelocities(std::vector<BoundaryPoint>&, MutableFloat&, const double, MersenneTwister&);
#endif

        //! Compute the modulus of the gradient of the signed distance function.
        void computeGradients();

        //! Calculate material area fraction enclosed by the discretised boundary.
        /*! \param boundary
                A reference to the discretised boundary.

            \return
                The total area fraction enclosed by the boundary.
         */
        double computeAreaFractions(const Boundary&);

        std::vector<double> signedDistance;     //!< The nodal signed distance function (level set).
        std::vector<double> velocity;           //!< The nodal normal velocity.
        std::vector<double> gradient;           //!< The nodal gradient of the level set function (modulus).
        std::vector<double> target;             //!< Signed distance target (for shape matching).
        std::vector<unsigned int> narrowBand;   //!< Indices of nodes in the narrow band.
        std::vector<unsigned int> mines;        //!< Indices of nodes at the edge of the narrow band.
        unsigned int nNarrowBand;               //!< The number of nodes in narrow band.
        unsigned int nMines;                    //!< The number of mine nodes.
        const double moveLimit;                 //!< The boundary movement limit (CFL condition).

        double area;                            //!< The total mesh area fraction enclosed by the boundary.
        Mesh mesh;                              //!< The fixed-grid mesh.

    private:
        unsigned int bandWidth;                 //!< The width of the narrow band region.
        bool isFixedDomain;                     //!< Whether the domain boundary is fixed.

        //! Default initialisation of the level set function (Swiss cheese configuration).
        void initialise();

        //! Initialise the level set from a vector of user-defined holes.
        /*! \param holes
                A vector of holes.
         */
        void initialise(const std::vector<Hole>&);

        //! Initialise the level set from a vector of user-defined points.
        /*! \param points
                A vector of point coordinates.
         */
        void initialise(const std::vector<Coord>&);

        //! Helper function for initialise methods.
        //! Initialises the level set function as the distance to the closest domain boundary.
        void closestDomainBoundary();

        //! Initialise the narrow band region.
        void initialiseNarrowBand();

        //! Initialise velocities for boundary nodes.
        /*! \param boundaryPoints
                A reference to a vector of boundary points.
         */
        void initialiseVelocities(const std::vector<BoundaryPoint>&);

        //! Compute the modulus of the gradient of the signed distance function at a node.
        /*! \param node
                The node index.

            \return
                The gradient at the node.
         */
        double computeGradient(const unsigned int) const;

        //! Compute Hamilton-Jacobi WENO gradient approximation.
        /*! \param v1
                The value of the function at the first stencil point.

            \param v2
                The value of the function at the second stencil point.

            \param v3
                The value of the function at the third stencil point.

            \param v4
                The value of the function at the fourth stencil point.

            \param v5
                The value of the function at the fifth stencil point.

            \return
                The smoothed function (gradient).
         */
        double gradHJWENO(double, double, double, double, double) const;

        //! Compute the minimum distance between a point and a line segment.
        /*! \param vertex1
                The coordinate of the first vertex.

            \param vertex2
                The coordinate of the second vertex.

            \param point
                The coordinates of the point of interest.

            \return
                The minimum distance to the line segment.
         */
        double pointToLineDistance(const Coord&, const Coord&, const Coord&) const;

        //! Test if a point lies inside a polygon.
        /*! \param point
                The coordinates of the point.

            \param vertices
                The vertices of the polygon (closed and ordered).

            \return
                Whether the point lies inside the polygon.
         */
        bool isInsidePolygon(const Coord&, const std::vector<Coord>&) const;

        //! Test if a point lies left, on, or right of an infinite line.
        /*! \param vertex1
                The coordinate of the first vertex.

            \param vertex2
                The coordinate of the second vertex.

            \param point
                The coordinates of the point of interest.

            \return
                > 0 if the point is left of the line through vertex1 and vertex2,
                = 0 if the point is on the line,
                < 0 if right of the line.
         */
        int isLeftOfLine(const Coord&, const Coord&, const Coord&) const;

        //! Calculate the material area for an element cut by the boundary.
        /*! \param element
                A reference to the element.

            \param boundary
                A reference to the discretised boundary.

            \return
                The area fraction.
         */
        double cutArea(const Element&, const Boundary&) const;

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
    };
}

#endif  /* _LEVELSET_H */
