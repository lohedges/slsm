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

#include "slsm.h"

int testBoundaryPoints()
{
    // Tests for the correct assignent of boundary points.
    //  1) Check that there are the correct number of boundary points and segments.
    //  2) Check that boundary point coordinates are correct.

    // Push hole into a vector container.
    std::vector<slsm::Hole> holes;
    holes.push_back(slsm::Hole(1, 1, 1));

    // Initialise a 1x1 level set domain.
    slsm::LevelSet levelSet(1, 1, holes);

    // Initialise the boundary object.
    slsm::Boundary boundary;

    // Set error number.
    errno = 0;

    // Sub test 1:
    // Horizontal boundary cutting through the middle of the element.

    // Place nodes on bottom edge of element outside zero contour.
    levelSet.signedDistance[0] = -1;
    levelSet.signedDistance[1] = -1;

    // Place nodes on top edge of element inside zero contour.
    levelSet.signedDistance[2] = 1;
    levelSet.signedDistance[3] = 1;

    // Discretise the boundary.
    boundary.discretise(levelSet);

    // Check the number of points and segments.
    slsm_check((boundary.nPoints == 2), "The number of boundary points is incorrect!");
    slsm_check((boundary.nSegments == 1), "The number of boundary segments is incorrect!");

    // Check the positions of the boundary points.

    // First point.
    slsm_check((std::abs(boundary.points[0].coord.x - 1) < 1e-6), "Position of boundary point is incorrect!");
    slsm_check((std::abs(boundary.points[0].coord.y - 0.5) < 1e-6), "Position of boundary point is incorrect!");

    // Second point.
    slsm_check((std::abs(boundary.points[1].coord.x) < 1e-6), "Position of boundary point is incorrect!");
    slsm_check((std::abs(boundary.points[1].coord.y - 0.5) < 1e-6), "Position of boundary point is incorrect!");

    // Sub test 2:
    // Vertical boundary cutting through the middle of the element.

    // Place nodes on left edge of element outside zero contour.
    levelSet.signedDistance[0] = -1;
    levelSet.signedDistance[2] = -1;

    // Place nodes on right edge of element inside zero contour.
    levelSet.signedDistance[1] = 1;
    levelSet.signedDistance[3] = 1;

    // Discretise the boundary.
    boundary.discretise(levelSet);

    // Check the number of points and segments.
    slsm_check((boundary.nPoints == 2), "The number of boundary points is incorrect!");
    slsm_check((boundary.nSegments == 1), "The number of boundary segments is incorrect!");

    // Check the positions of the boundary points.

    // First point.
    slsm_check((std::abs(boundary.points[0].coord.x - 0.5) < 1e-6), "Position of boundary point is incorrect!");
    slsm_check((std::abs(boundary.points[0].coord.y) < 1e-6), "Position of boundary point is incorrect!");

    // Second point.
    slsm_check((std::abs(boundary.points[1].coord.x - 0.5) < 1e-6), "Position of boundary point is incorrect!");
    slsm_check((std::abs(boundary.points[1].coord.y - 1) < 1e-6), "Position of boundary point is incorrect!");

    // Sub test 3:
    // Two boundary segements cutting through upper and lower triangles.

    // Place bottom left and top right nodes outside zero contour.
    levelSet.signedDistance[0] = -1;
    levelSet.signedDistance[3] = -1;

    // Place bottom right and top left nodes inside zero contour.
    levelSet.signedDistance[1] = 1;
    levelSet.signedDistance[2] = 1;

    // Discretise the boundary.
    boundary.discretise(levelSet);

    // Check the number of points and segments.
    slsm_check((boundary.nPoints == 4), "The number of boundary points is incorrect!");
    slsm_check((boundary.nSegments == 2), "The number of boundary segments is incorrect!");

    // Check the positions of the boundary points.

    // First point.
    slsm_check((std::abs(boundary.points[0].coord.x - 0.5) < 1e-6), "Position of boundary point is incorrect!");
    slsm_check((std::abs(boundary.points[0].coord.y) < 1e-6), "Position of boundary point is incorrect!");

    // Second point.
    slsm_check((std::abs(boundary.points[1].coord.x - 1) < 1e-6), "Position of boundary point is incorrect!");
    slsm_check((std::abs(boundary.points[1].coord.y - 0.5) < 1e-6), "Position of boundary point is incorrect!");

    // Third point.
    slsm_check((std::abs(boundary.points[2].coord.x - 0.5) < 1e-6), "Position of boundary point is incorrect!");
    slsm_check((std::abs(boundary.points[2].coord.y - 1) < 1e-6), "Position of boundary point is incorrect!");

    // Fourth point.
    slsm_check((std::abs(boundary.points[3].coord.x) < 1e-6), "Position of boundary point is incorrect!");
    slsm_check((std::abs(boundary.points[3].coord.y - 0.5) < 1e-6), "Position of boundary point is incorrect!");

    // Sub test 4:
    // An element with a single segment formed from the diagonal between
    // the bottom left and top right nodes (nodes 0 and 3).

    // Place bottom left and top right nodes on zero contour.
    levelSet.signedDistance[0] = 0;
    levelSet.signedDistance[3] = 0;

    // Place bottom right node outside.
    levelSet.signedDistance[1] = -1;

    // Place top left node inside.
    levelSet.signedDistance[2] = 1;

    // Discretise the boundary.
    boundary.discretise(levelSet);

    // Check the positions of the boundary points.

    // First point.
    slsm_check((std::abs(boundary.points[0].coord.x - levelSet.mesh.nodes[0].coord.x) < 1e-6), "Position of boundary point is incorrect!");
    slsm_check((std::abs(boundary.points[0].coord.y - levelSet.mesh.nodes[0].coord.y) < 1e-6), "Position of boundary point is incorrect!");

    // Second point.
    slsm_check((std::abs(boundary.points[1].coord.x - levelSet.mesh.nodes[3].coord.x) < 1e-6), "Position of boundary point is incorrect!");
    slsm_check((std::abs(boundary.points[1].coord.y - levelSet.mesh.nodes[3].coord.y) < 1e-6), "Position of boundary point is incorrect!");

    return 0;

error:
    return 1;
}

int testBoundarySegments()
{
    // Tests for the correct assignent of boundary segments.
    //  1) Check for correct indexing of start and end points.
    //  2) Check that segments are assigned to the correct element.

    // Push hole into a vector container.
    std::vector<slsm::Hole> holes;
    holes.push_back(slsm::Hole(1, 1, 1));

    // Initialise a 1x1 level set domain.
    slsm::LevelSet levelSet(1, 1, holes);

    // Initialise the boundary object.
    slsm::Boundary boundary;

    // Set error number.
    errno = 0;

    // Sub test 1:
    // Boundary cutting horizontally through the middle of the element.

    // Place nodes on bottom edge of element outside zero contour.
    levelSet.signedDistance[0] = -1;
    levelSet.signedDistance[1] = -1;

    // Place nodes on top edge of element inside zero contour.
    levelSet.signedDistance[2] = 1;
    levelSet.signedDistance[3] = 1;

    // Discretise the boundary.
    boundary.discretise(levelSet);

    // Check the number of points and segments.
    slsm_check((boundary.nPoints == 2), "The number of boundary points is incorrect!");
    slsm_check((boundary.nSegments == 1), "The number of boundary segments is incorrect!");

    // Check the boundary segment.

    // Start point.
    slsm_check((boundary.segments[0].start == 0), "Start point of segment 0 is incorrect!");

    // End point.
    slsm_check((boundary.segments[0].end == 1), "End point of segment 0 is incorrect!");

    // Element.
    slsm_check((boundary.segments[0].element == 0), "Element of segment 0 is incorrect!");

    // Sub test 2:
    // Boundary segment running from zero contour at lower left node to a boundary
    // point at the middle of the right-hand element edge.

    // Place bottom left node on zero contour.
    levelSet.signedDistance[0] = 0;

    // Place bottom right node on outside zero contour.
    levelSet.signedDistance[1] = -1;

    // Place nodes on top edge of element inside zero contour.
    levelSet.signedDistance[2] = 1;
    levelSet.signedDistance[3] = 1;

    // Discretise the boundary.
    boundary.discretise(levelSet);

    // Check the number of points and segments.
    slsm_check((boundary.nPoints == 2), "The number of boundary points is incorrect!");
    slsm_check((boundary.nSegments == 1), "The number of boundary segments is incorrect!");

    // Check the boundary segment.
    // The start point will be the boundary point lying between nodes 1 and 2 on
    // the right edge of the element. The end point of the segment should lie on
    // the lower left node of the grid.

    // Start point.
    slsm_check((boundary.segments[0].start == 0), "Start point of segment 0 is incorrect!");

    // End point.
    slsm_check((boundary.segments[0].end == 1), "End point of segment 0 is incorrect!");

    // Element.
    slsm_check((boundary.segments[0].element == 0), "Element of segment 0 is incorrect!");

    // Sub test 2:
    // An element with a single segment formed from the diagonal between
    // the bottom left and top right nodes (nodes 0 and 3).

    // Place bottom left and top right nodes on zero contour.
    levelSet.signedDistance[0] = 0;
    levelSet.signedDistance[3] = 0;

    // Place bottom right node outside.
    levelSet.signedDistance[1] = -1;

    // Place top left node inside.
    levelSet.signedDistance[2] = 1;

    // Discretise the boundary.
    boundary.discretise(levelSet);

    // Check the boundary segment. The single boundary segment should lie
    // along the diagonal between nodes 0 and 3 (which lie on the zero contour).

    slsm_check((boundary.nPoints == 2), "The number of boundary points is incorrect!");
    slsm_check((boundary.nSegments == 1), "The number of boundary segments is incorrect!");

    // Start point.
    slsm_check((boundary.segments[0].start == 0), "Start point of segment 0 is incorrect!");

    // End point.
    slsm_check((boundary.segments[0].end == 1), "End point of segment 0 is incorrect!");

    // Element.
    slsm_check((boundary.segments[0].element == 0), "Element of segment 0 is incorrect!");

    return 0;

error:
    return 1;
}

int testBoundarySymmetry()
{
    // Push hole into a vector container.
    std::vector<slsm::Hole> holes;
    holes.push_back(slsm::Hole(1, 1, 1));

    // Initialise a 1x1 level set domain.
    slsm::LevelSet levelSet(1, 1, holes);

    // Initialise the boundary object.
    slsm::Boundary boundary;

    // Set error number.
    errno = 0;

    // Place bottom left and top right nodes outside zero contour.
    levelSet.signedDistance[0] = -1;
    levelSet.signedDistance[3] = -1;

    // Place bottom right and top left nodes inside zero contour.
    levelSet.signedDistance[1] = 1;
    levelSet.signedDistance[2] = 1;

    // Discretise the boundary.
    boundary.discretise(levelSet);

    // Boundary should cut diagonally through the upper and lower triangles of the element.

    // Check the number of points and segments.
    slsm_check((boundary.nPoints == 4), "The number of boundary points is incorrect!");
    slsm_check((boundary.nSegments == 2), "The number of boundary segments is incorrect!");

    // Check the positions of the boundary points.

    // First point.
    slsm_check((std::abs(boundary.points[0].coord.x - 0.5) < 1e-6), "Position of boundary point is incorrect!");
    slsm_check((std::abs(boundary.points[0].coord.y) < 1e-6), "Position of boundary point is incorrect!");

    // Second point.
    slsm_check((std::abs(boundary.points[1].coord.x - 1) < 1e-6), "Position of boundary point is incorrect!");
    slsm_check((std::abs(boundary.points[1].coord.y - 0.5) < 1e-6), "Position of boundary point is incorrect!");

    // Third point.
    slsm_check((std::abs(boundary.points[2].coord.x - 0.5) < 1e-6), "Position of boundary point is incorrect!");
    slsm_check((std::abs(boundary.points[2].coord.y - 1) < 1e-6), "Position of boundary point is incorrect!");

    // Fourth point.
    slsm_check((std::abs(boundary.points[3].coord.x) < 1e-6), "Position of boundary point is incorrect!");
    slsm_check((std::abs(boundary.points[3].coord.y - 0.5) < 1e-6), "Position of boundary point is incorrect!");

    // Invert the signed distance function.
    for (unsigned int i=0;i<4;i++)
        levelSet.signedDistance[i] *= -1;

    // Discretise the boundary.
    boundary.discretise(levelSet);

    // Boundary points should be unchanged.

    // Check the number of points and segments.
    slsm_check((boundary.nPoints == 4), "The number of boundary points is incorrect!");
    slsm_check((boundary.nSegments == 2), "The number of boundary segments is incorrect!");

    // First point.
    slsm_check((std::abs(boundary.points[0].coord.x - 0.5) < 1e-6), "Position of boundary point is incorrect!");
    slsm_check((std::abs(boundary.points[0].coord.y) < 1e-6), "Position of boundary point is incorrect!");

    // Second point.
    slsm_check((std::abs(boundary.points[1].coord.x - 1) < 1e-6), "Position of boundary point is incorrect!");
    slsm_check((std::abs(boundary.points[1].coord.y - 0.5) < 1e-6), "Position of boundary point is incorrect!");

    // Third point.
    slsm_check((std::abs(boundary.points[2].coord.x - 0.5) < 1e-6), "Position of boundary point is incorrect!");
    slsm_check((std::abs(boundary.points[2].coord.y - 1) < 1e-6), "Position of boundary point is incorrect!");

    // Fourth point.
    slsm_check((std::abs(boundary.points[3].coord.x) < 1e-6), "Position of boundary point is incorrect!");
    slsm_check((std::abs(boundary.points[3].coord.y - 0.5) < 1e-6), "Position of boundary point is incorrect!");

    return 0;

error:
    return 1;
}

int testConnectivity()
{
    // Push hole into a vector container.
    std::vector<slsm::Hole> holes;
    holes.push_back(slsm::Hole(1, 1, 1));

    // Initialise a 1x1 level set domain.
    slsm::LevelSet levelSet(1, 1, holes);

    // Initialise the boundary object.
    slsm::Boundary boundary;

    // Set error number.
    errno = 0;

    // Place bottom left and top right nodes outside zero contour.
    levelSet.signedDistance[0] = -1;
    levelSet.signedDistance[3] = -1;

    // Place bottom right and top left nodes inside zero contour.
    levelSet.signedDistance[1] = 1;
    levelSet.signedDistance[2] = 1;

    // Discretise the boundary.
    boundary.discretise(levelSet);

    // Boundary should cut diagonally through the upper and lower triangles of the element.

    // Check connectivity between nodes and boundary points.

    // First node.
    slsm_check((levelSet.mesh.nodes[0].nBoundaryPoints == 2), "Incorrect number of boundary points associated with node 0!");
    slsm_check((levelSet.mesh.nodes[0].boundaryPoints[0] == 0), "Boundary point 0 of node 0 is incorrect!");
    slsm_check((levelSet.mesh.nodes[0].boundaryPoints[1] == 3), "Boundary point 1 of node 0 is incorrect!");

    // Second node.
    slsm_check((levelSet.mesh.nodes[1].nBoundaryPoints == 2), "Incorrect number of boundary points associated with node 1!");
    slsm_check((levelSet.mesh.nodes[1].boundaryPoints[0] == 0), "Boundary point 0 of node 1 is incorrect!");
    slsm_check((levelSet.mesh.nodes[1].boundaryPoints[1] == 1), "Boundary point 1 of node 1 is incorrect!");

    // Third node.
    slsm_check((levelSet.mesh.nodes[2].nBoundaryPoints == 2), "Incorrect number of boundary points associated with node 2!");
    slsm_check((levelSet.mesh.nodes[2].boundaryPoints[0] == 2), "Boundary point 0 of node 2 is incorrect!");
    slsm_check((levelSet.mesh.nodes[2].boundaryPoints[1] == 3), "Boundary point 1 of node 2 is incorrect!");

    // Fourth node.
    slsm_check((levelSet.mesh.nodes[3].nBoundaryPoints == 2), "Incorrect number of boundary points associated with node 3!");
    slsm_check((levelSet.mesh.nodes[3].boundaryPoints[0] == 1), "Boundary point 0 of node 3 is incorrect!");
    slsm_check((levelSet.mesh.nodes[3].boundaryPoints[1] == 2), "Boundary point 1 of node 3 is incorrect!");

    // Check boundary segment connectivity.
    slsm_check((levelSet.mesh.elements[0].nBoundarySegments == 2), "Incorrect number of boundary segments for node 0!");
    slsm_check((levelSet.mesh.elements[0].boundarySegments[0] == 0), "Boundary segment 0 of element 0 is incorrect!");
    slsm_check((levelSet.mesh.elements[0].boundarySegments[1] == 1), "Boundary segment 1 of element 0 is incorrect!");

    return 0;

error:
    return 1;
}

int testAreaFraction()
{
    // Push hole into a vector container.
    std::vector<slsm::Hole> holes;
    holes.push_back(slsm::Hole(1, 1, 1));

    // Initialise a 1x1 level set domain.
    slsm::LevelSet levelSet(1, 1, holes);

    // Initialise the boundary object.
    slsm::Boundary boundary;

    // Set error number.
    errno = 0;

    // Sub test 1:
    // Horizontal boundary cutting through the middle of the element.

    // Place nodes on bottom edge of element outside zero contour.
    levelSet.signedDistance[0] = -1;
    levelSet.signedDistance[1] = -1;

    // Place nodes on top edge of element inside zero contour.
    levelSet.signedDistance[2] = 1;
    levelSet.signedDistance[3] = 1;

    // Discretise the boundary.
    boundary.discretise(levelSet);

    // Calculate the element area.
    levelSet.computeAreaFractions(boundary);

    // Check that the element is half filled.
    slsm_check((levelSet.mesh.elements[0].area == 0.5), "Element area fraction is incorrect!");

    // Sub test 2:
    // Vertical boundary cutting through the middle of the element.

    // Place nodes on left edge of element outside zero contour.
    levelSet.signedDistance[0] = -1;
    levelSet.signedDistance[2] = -1;

    // Place nodes on right edge of element inside zero contour.
    levelSet.signedDistance[1] = 1;
    levelSet.signedDistance[3] = 1;

    // Discretise the boundary.
    boundary.discretise(levelSet);

    // Calculate the element area.
    levelSet.computeAreaFractions(boundary);

    // Check that the element is half filled.
    slsm_check((levelSet.mesh.elements[0].area == 0.5), "Element area fraction is incorrect!");

    // Sub test 3:
    // Zero contour along diagonal.

    // Place bottom left and top right nodes on zero contour.
    levelSet.signedDistance[0] = 0;
    levelSet.signedDistance[3] = 0;

    // Place bottom right node inside and top left node outside.
    levelSet.signedDistance[1] = 1;
    levelSet.signedDistance[2] = -1;

    // Discretise the boundary.
    boundary.discretise(levelSet);

    // Calculate the element area.
    levelSet.computeAreaFractions(boundary);

    // Check that the element is half filled.
    slsm_check((levelSet.mesh.elements[0].area == 0.5), "Element area fraction is incorrect!");

    // Sub test 4:
    // Zero contour along diagonal in top right quadrant.

    // Place bottom left, bottom right, and top left nodes inside.
    levelSet.signedDistance[0] = 1;
    levelSet.signedDistance[1] = 1;
    levelSet.signedDistance[2] = 1;

    // Place top right node outside.
    levelSet.signedDistance[3] = -1;

    // Discretise the boundary.
    boundary.discretise(levelSet);

    // Calculate the element area.
    levelSet.computeAreaFractions(boundary);

    // Check the element area.
    // Full element minus half a quarter element, i.e. 1.0 - 1/8 = 0.875.
    slsm_check((levelSet.mesh.elements[0].area == 0.875), "Element area fraction is incorrect!");

    // Sub test 4:
    // Invert the level set variables from sub-test 3.
    // The element area should now be one minus the answer above, i.e. 1/8.

    // Place bottom left, bottom right, and top left nodes outside.
    levelSet.signedDistance[0] = -1;
    levelSet.signedDistance[1] = -1;
    levelSet.signedDistance[2] = -1;

    // Place top right node inside.
    levelSet.signedDistance[3] = 1;

    // Discretise the boundary.
    boundary.discretise(levelSet);

    // Calculate the element area.
    levelSet.computeAreaFractions(boundary);

    // Check the element area.
    slsm_check((levelSet.mesh.elements[0].area == (1.0/8.0)), "Element area fraction is incorrect!");

    return 0;

error:
    return 1;
}

int all_tests()
{
    mu_suite_start();

    mu_run_test(testBoundaryPoints);
    mu_run_test(testBoundarySegments);
    mu_run_test(testBoundarySymmetry);
    mu_run_test(testConnectivity);
    mu_run_test(testAreaFraction);

    return 0;
}

RUN_TESTS(all_tests);
