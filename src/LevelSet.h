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

#ifndef _LEVELSET_H
#define _LEVELSET_H

#include <cmath>
#include <cstdlib>

#include "Common.h"
#include "FastMarchingMethod.h"
#include "Hole.h"
#include "MersenneTwister.h"
#include "Mesh.h"

/*! \file LevelSet.h
    \brief A class for the level set function.
 */

namespace lsm
{
    //! Forward declaration of BoundaryPoint data structure.
    //! Might want to move all non-class types to Common.h, i.e. affiliated types.
    struct BoundaryPoint;

    /*! \brief A class for the level set function.

        The level set is represented as a signed distance function from the zero
        countour. Positive values are inside the structure, negative values are
        outside.

        The class provides methods for initialising and reinitialising the signed
        distance function. The default initialisation uses a set of linearly spaced
        circular holes to create the classic "Swiss cheese" starting configuration.
        Alternatively, the user can pass the constructor a vector of holes that
        will be used to define the initial structure.

        Reinitialisation of the signed distance function is performed using an
        implementation of the Fast Marching Method. A second-order stencil is used
        for the upwind finite difference scheme where possible.

        Functionality is also provided for tracking nodes that are part of the
        narrow band region around the zero contour, as well as mine nodes at
        the edge of the narrow band.

        Note that the level set grid need not be the same resolution as the
        finite-element mesh. However, for the simple two-dimensional problems
        considered here we will use a grid of the same size.
     */
    class LevelSet
    {
    public:
        //! Constructor.
        /*! \param mesh_
                A reference to the fixed grid, finite-element mesh.

            \param moveLimit_
                The CFL limit (in units of the mesh grid spacing).

            \param bandWidth_
                The width of the narrow band region.

            \param isFixed_
                Whether the domain boundary is fixed.
         */
        LevelSet(Mesh&, double moveLimit_ = 0.5, unsigned int bandWidth_ = 6, bool isFixed_ = false);

        //! Constructor.
        /*! \param mesh_
                A reference to the fixed grid, finite-element mesh.

            \param holes
                A vector of holes.

            \param moveLimit_
                The CFL limit (in units of the mesh grid spacing).

            \param bandWidth_
                The width of the narrow band region.

            \param isFixed_
                Whether the domain boundary is fixed.
         */
        LevelSet(Mesh&, const std::vector<Hole>&, double moveLimit_ = 0.5,
            unsigned int bandWidth_ = 6, bool isFixed_ = false);

        //! Update the level set function.
        /*! \param timeStep
                The time step.

            \return
                Whether the signed distance was reinitialised.
         */
        bool update(double);

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

        //! Compute the modulus of the gradient of the signed distance function.
        void computeGradients();

        std::vector<double> signedDistance;     //!< The nodal signed distance function (level set).
        std::vector<double> velocity;           //!< The nodal signed normal velocity.
        std::vector<double> gradient;           //!< The nodal gradient of the level set function (modulus).
        std::vector<unsigned int> narrowBand;   //!< Indices of nodes in the narrow band.
        std::vector<unsigned int> mines;        //!< Indices of nodes at the edge of the narrow band.
        const unsigned int nNodes;              //!< The number of nodes in the finite element grid.
        unsigned int nNarrowBand;               //!< The number of nodes in narrow band.
        unsigned int nMines;                    //!< The number of mine nodes.
        const double moveLimit;                 //!< The boundary movement limit (CFL condition).

    private:
        Mesh& mesh;                             //!< A reference to the finite element mesh.
        unsigned int bandWidth;                 //!< The width of the narrow band region.
        bool isFixed;                           //!< Whether the domain boundary is fixed.

        //! Default initialisation of the level set function (Swiss cheese configuration).
        void initialise();

        //! Initialise the level set from a vector of user-defined holes.
        /* \param holes
                A vector of holes.
         */
        void initialise(const std::vector<Hole>&);

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
        double computeGradient(const unsigned int);

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
        double gradHJWENO(double, double, double, double, double);
    };
}

#endif  /* _LEVELSET_H */
