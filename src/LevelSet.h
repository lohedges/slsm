/*
 * File:	LevelSet.h
 * Author:	lester
 */

#ifndef _LEVELSET_H
#define _LEVELSET_H

#include <algorithm>

#include "Common.h"
#include "Hole.h"
#include "Mesh.h"
#include "FastMarchingMethod.h"

/*! \file LevelSet.h
    \brief A class for the level set function.
*/

class LevelSet
{
public:
    //! Constructor.
    /*! \param mesh_
            A reference to the fixed grid, finite-element mesh.

        \param bandWidth_
            The width of the narrow band region.
     */
	LevelSet(Mesh&, unsigned int);

    //! Constructor.
    /*! \param mesh_
            A reference to the fixed grid, finite-element mesh.

        \param bandWidth_
            The width of the narrow band region.

        \param holes
            A vector of holes.
     */
	LevelSet(Mesh&, unsigned int, const std::vector<Hole>&);

    //! Update the level set function.
    void update();

    //! Re-initialise the level set function to a signed distance function using
    //! the fast marching method.
    void reinitialise();

    std::vector<double> signedDistance;     //!< The nodal signed distance function (level set).
    std::vector<double> gradient;           //!< The nodal gradient of the level set function.
    std::vector<unsigned int> narrowBand;   //!< Indices of nodes in the narrow band.
    std::vector<unsigned int> mines;        //!< Indices of nodes at the edge of the narrow band.
    const unsigned int nNodes;              //!< The number of nodes in the finite element grid.
    unsigned int nNarrowBand;               //!< The number of nodes in narrow band.
    unsigned int nMines;                    //!< The number of mine nodes.

private:
    Mesh& mesh;                             //!< A reference to the finite element mesh.
    unsigned int bandWidth;                 //!< The width of the narrow band region.
    unsigned int boundaryNode;              //!< The node that is currently closest to the zero isocontour.
    double minDistance;                     //!< Absolute minimum value of the signed distance function.

    //! Default initialisation of the level set function (Swiss cheese configuration).
    void initialise();

    //! Initialise the level set from a vector of user-defined holes.
    /*/
        \param holes
            A vector of holes.
     */
    void initialise(const std::vector<Hole>&);

    //! Helper function for initialise methods.
    //! Initialises the level set function as the distance to the closest domain boundary.
    void closestDomainBoundary();

    //! Initialise the narrow band region.
    void initialiseNarrowBand();
};

#endif	/* _LEVELSET_H */
