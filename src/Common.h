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

#ifndef _COMMON_H
#define _COMMON_H

/*! \file Common.h
    \brief Common data types.
 */

namespace slsm
{
    //! Two-dimensional coordinate information.
    struct Coord
    {
        //! Constructor.
        Coord() {};

        //! Constructor.
        /*! \param x_
                The x coordinate.

            \param y_
                The y coordinate.
         */
        Coord(double x_, double y_) {x = x_; y = y_;};

        double x;   //!< The x coordinate.
        double y;   //!< The y coordinate.
    };
}

#endif  /* _COMMON_H */
