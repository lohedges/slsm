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

#include "Hole.h"

/*! \file Hole.cpp
    \brief A simple circular hole data type.
 */

namespace slsm
{
    Hole::Hole() {}

    Hole::Hole(double x, double y, double r) : coord(x, y), r(r)
    {
    }

    Hole::Hole(Coord& coord_, double r) : coord(coord_), r(r)
    {
    }
}
