/*
 * File:	Hole.h
 * Author:	lester
 */

#ifndef _HOLE_H
#define _HOLE_H

#include "Common.h"

//! Class for handling circular holes.
class Hole
{
public:
    //! Default constructor.
    Hole();

    //! Constructor.
    /*! \param x
            The x coordinate of the hole.

        \param y
            The x coordinate of the hole.

        \param r
            The radius of the hole.
     */
    Hole(double, double, double);

    //! Constructor.
    /*! \param coord_
            The x-y coordinates of the hole.

        \param r
            The radius of the hole.
     */
    Hole(Coord&, double);

    Coord coord;                            //!< Coordinates.
    double r;                               //!< Radius.
};

#endif	/* _HOLE_H */
