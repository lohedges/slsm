/*
 * File:    Hole.cpp
 * Author:  lester
 */

#include "Hole.h"

Hole::Hole() {}

Hole::Hole(double x, double y, double r) : r(r)
{
    coord.x = x;
    coord.y = y;
}

Hole::Hole(Coord& coord_, double r) : r(r)
{
    coord.x = coord_.x;
    coord.y = coord_.y;
}
