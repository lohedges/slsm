/*
 * File:	Boundary.h
 * Author:	lester
 */

#ifndef _BOUNDARY_H
#define _BOUNDARY_H

struct BoundarySegment
{
    unsigned int node1;     //!< Index for first node end point.
    unsigned int node2;     //!< Index for second node end point.
    unsigned int element;   //!< The element cut by the boundary segment.
    double length;          //!< Length of the boundary segment.
    double weight;          //!< Weighting factor for boundary segment.
};

class Boundary
{
public:
	Boundary();
	~Boundary();
};

#endif	/* _BOUNDARY_H */
