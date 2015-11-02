/*
 * File:   boundary_tests.cpp
 * Author: lester
 */

#include "lsm.h"

int testBoundaryPoints()
{
    // Initialise a 1x1 non-periodic mesh.
    Mesh mesh(1, 1, false);

    // Push hole into a vector container.
    std::vector<Hole> holes;
    holes.push_back(Hole(1, 1, 1));

    // Initialise the level set object.
    LevelSet levelSet(mesh, 3, holes);

    // Initialise the boundary object.
    Boundary boundary(mesh, levelSet);

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
    boundary.discretise();

    // Check the number of points and segments.
    check((boundary.nPoints == 2), "The number of boundary points is incorrect!");
    check((boundary.nSegments == 1), "The number of boundary segments is incorrect!");

    // Check the positions of the boundary points.

    // First point.
    check((std::abs(boundary.points[0].x - 1) < 1e-6), "Position of boundary point is incorrect!");
    check((std::abs(boundary.points[0].y - 0.5) < 1e-6), "Position of boundary point is incorrect!");

    // Second point.
    check((std::abs(boundary.points[1].x) < 1e-6), "Position of boundary point is incorrect!");
    check((std::abs(boundary.points[1].y - 0.5) < 1e-6), "Position of boundary point is incorrect!");

    // Sub test 2:
    // Vertical boundary cutting through the middle of the element.

    // Place nodes on left edge of element outside zero contour.
    levelSet.signedDistance[0] = -1;
    levelSet.signedDistance[2] = -1;

    // Place nodes on right edge of element inside zero contour.
    levelSet.signedDistance[1] = 1;
    levelSet.signedDistance[3] = 1;

    // Discretise the boundary.
    boundary.discretise();

    // Check the number of points and segments.
    check((boundary.nPoints == 2), "The number of boundary points is incorrect!");
    check((boundary.nSegments == 1), "The number of boundary segments is incorrect!");

    // Check the positions of the boundary points.

    // First point.
    check((std::abs(boundary.points[0].x - 0.5) < 1e-6), "Position of boundary point is incorrect!");
    check((std::abs(boundary.points[0].y) < 1e-6), "Position of boundary point is incorrect!");

    // Second point.
    check((std::abs(boundary.points[1].x - 0.5) < 1e-6), "Position of boundary point is incorrect!");
    check((std::abs(boundary.points[1].y - 1) < 1e-6), "Position of boundary point is incorrect!");

    // Sub test 3:
    // Two boundary segements cutting through upper and lower triangles.

    // Place bottom left and top right nodes outside zero contour.
    levelSet.signedDistance[0] = -1;
    levelSet.signedDistance[3] = -1;

    // Place bottom right and top left nodes inside zero contour.
    levelSet.signedDistance[1] = 1;
    levelSet.signedDistance[2] = 1;

    // Discretise the boundary.
    boundary.discretise();

    // Check the number of points and segments.
    check((boundary.nPoints == 4), "The number of boundary points is incorrect!");
    check((boundary.nSegments == 2), "The number of boundary segments is incorrect!");

    // Check the positions of the boundary points.

    // First point.
    check((std::abs(boundary.points[0].x - 0.5) < 1e-6), "Position of boundary point is incorrect!");
    check((std::abs(boundary.points[0].y) < 1e-6), "Position of boundary point is incorrect!");

    // Second point.
    check((std::abs(boundary.points[1].x - 1) < 1e-6), "Position of boundary point is incorrect!");
    check((std::abs(boundary.points[1].y - 0.5) < 1e-6), "Position of boundary point is incorrect!");

    // Third point.
    check((std::abs(boundary.points[2].x - 0.5) < 1e-6), "Position of boundary point is incorrect!");
    check((std::abs(boundary.points[2].y - 1) < 1e-6), "Position of boundary point is incorrect!");

    // Fourth point.
    check((std::abs(boundary.points[3].x) < 1e-6), "Position of boundary point is incorrect!");
    check((std::abs(boundary.points[3].y - 0.5) < 1e-6), "Position of boundary point is incorrect!");

    return 0;

error:
    return 1;
}

int testBoundarySegments()
{
    // Initialise a 1x1 non-periodic mesh.
    Mesh mesh(1, 1, false);

    // Push hole into a vector container.
    std::vector<Hole> holes;
    holes.push_back(Hole(1, 1, 1));

    // Initialise the level set object.
    LevelSet levelSet(mesh, 3, holes);

    // Initialise the boundary object.
    Boundary boundary(mesh, levelSet);

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
    boundary.discretise();

    // Check the number of points and segments.
    check((boundary.nPoints == 2), "The number of boundary points is incorrect!");
    check((boundary.nSegments == 1), "The number of boundary segments is incorrect!");

    // Check the boundary segment.

    // N.B.
    // All start and end points should be boundary points, so the indexing is
    // shifted by nNodes.

    // Start point.
    check((boundary.segments[0].start == mesh.nNodes), "Start point of segment 0 is incorrect!");

    // End point.
    check((boundary.segments[0].end == (mesh.nNodes + 1)), "End point of segment 0 is incorrect!");

    // Element.
    check((boundary.segments[0].element == 0), "Element of segment 0 is incorrect!");

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
    boundary.discretise();

    // Check the number of points and segments.
    check((boundary.nPoints == 1), "The number of boundary points is incorrect!");
    check((boundary.nSegments == 1), "The number of boundary segments is incorrect!");

    // Check the boundary segment.
    // The start point will be the boundary point lying between nodes 1 and 2 on
    // the right edge of the element.
    // The end point of the segment should be the lower left node of the grid.

    // Start point.
    check((boundary.segments[0].start == mesh.nNodes), "Start point of segment 0 is incorrect!");

    // End point.
    check((boundary.segments[0].end == 0), "End point of segment 0 is incorrect!");

    // Element.
    check((boundary.segments[0].element == 0), "Element of segment 0 is incorrect!");

    // Sub test 2:
    // An element with no boundary points and a single segment formed from
    // the diagonal between the bottom left and top right nodes.

    // Place bottom left and top right nodes on zero contour.
    levelSet.signedDistance[0] = 0;
    levelSet.signedDistance[3] = 0;

    // Place bottom right node outside.
    levelSet.signedDistance[2] = -1;

    // Place top left node inside.
    levelSet.signedDistance[2] = 1;

    // Discretise the boundary.
    boundary.discretise();

    // Check the boundary segment.
    // There should be no cut edges, i.e. no boundary points.
    // The single boundary segment should lie along the diagonal between
    // nodes 0 and 3 (which lie on the zero contour).

    check((boundary.nPoints == 0), "The number of boundary points is incorrect!");
    check((boundary.nSegments == 1), "The number of boundary segments is incorrect!");

    // Start point.
    check((boundary.segments[0].start == 0), "Start point of segment 0 is incorrect!");

    // End point.
    check((boundary.segments[0].end == (mesh.nNodes - 1)), "End point of segment 0 is incorrect!");

    // Element.
    check((boundary.segments[0].element == 0), "Element of segment 0 is incorrect!");

    return 0;

error:
    return 1;
}

int testBoundarySymmetry()
{
    // Initialise a 1x1 non-periodic mesh.
    Mesh mesh(1, 1, false);

    // Push hole into a vector container.
    std::vector<Hole> holes;
    holes.push_back(Hole(1, 1, 1));

    // Initialise the level set object.
    LevelSet levelSet(mesh, 3, holes);

    // Initialise the boundary object.
    Boundary boundary(mesh, levelSet);

    // Set error number.
    errno = 0;

    // Place bottom left and top right nodes outside zero contour.
    levelSet.signedDistance[0] = -1;
    levelSet.signedDistance[3] = -1;

    // Place bottom right and top left nodes inside zero contour.
    levelSet.signedDistance[1] = 1;
    levelSet.signedDistance[2] = 1;

    // Discretise the boundary.
    boundary.discretise();

    // Boundary should cut diagonally through the upper and lower triangles of the element.

    // Check the number of points and segments.
    check((boundary.nPoints == 4), "The number of boundary points is incorrect!");
    check((boundary.nSegments == 2), "The number of boundary segments is incorrect!");

    // Check the positions of the boundary points.

    // First point.
    check((std::abs(boundary.points[0].x - 0.5) < 1e-6), "Position of boundary point is incorrect!");
    check((std::abs(boundary.points[0].y) < 1e-6), "Position of boundary point is incorrect!");

    // Second point.
    check((std::abs(boundary.points[1].x - 1) < 1e-6), "Position of boundary point is incorrect!");
    check((std::abs(boundary.points[1].y - 0.5) < 1e-6), "Position of boundary point is incorrect!");

    // Third point.
    check((std::abs(boundary.points[2].x - 0.5) < 1e-6), "Position of boundary point is incorrect!");
    check((std::abs(boundary.points[2].y - 1) < 1e-6), "Position of boundary point is incorrect!");

    // Fourth point.
    check((std::abs(boundary.points[3].x) < 1e-6), "Position of boundary point is incorrect!");
    check((std::abs(boundary.points[3].y - 0.5) < 1e-6), "Position of boundary point is incorrect!");

    // Invert the signed distance function.
    for (unsigned int i=0;i<4;i++)
        levelSet.signedDistance[i] *= -1;

    // Discretise the boundary.
    boundary.discretise();

    // Boundary points should be unchanged.

    // Check the number of points and segments.
    check((boundary.nPoints == 4), "The number of boundary points is incorrect!");
    check((boundary.nSegments == 2), "The number of boundary segments is incorrect!");

    // First point.
    check((std::abs(boundary.points[0].x - 0.5) < 1e-6), "Position of boundary point is incorrect!");
    check((std::abs(boundary.points[0].y) < 1e-6), "Position of boundary point is incorrect!");

    // Second point.
    check((std::abs(boundary.points[1].x - 1) < 1e-6), "Position of boundary point is incorrect!");
    check((std::abs(boundary.points[1].y - 0.5) < 1e-6), "Position of boundary point is incorrect!");

    // Third point.
    check((std::abs(boundary.points[2].x - 0.5) < 1e-6), "Position of boundary point is incorrect!");
    check((std::abs(boundary.points[2].y - 1) < 1e-6), "Position of boundary point is incorrect!");

    // Fourth point.
    check((std::abs(boundary.points[3].x) < 1e-6), "Position of boundary point is incorrect!");
    check((std::abs(boundary.points[3].y - 0.5) < 1e-6), "Position of boundary point is incorrect!");

    return 0;

error:
    return 1;
}

int testConnectivity()
{
    // Initialise a 1x1 non-periodic mesh.
    Mesh mesh(1, 1, false);

    // Push hole into a vector container.
    std::vector<Hole> holes;
    holes.push_back(Hole(1, 1, 1));

    // Initialise the level set object.
    LevelSet levelSet(mesh, 3, holes);

    // Initialise the boundary object.
    Boundary boundary(mesh, levelSet);

    // Set error number.
    errno = 0;

    // Place bottom left and top right nodes outside zero contour.
    levelSet.signedDistance[0] = -1;
    levelSet.signedDistance[3] = -1;

    // Place bottom right and top left nodes inside zero contour.
    levelSet.signedDistance[1] = 1;
    levelSet.signedDistance[2] = 1;

    // Discretise the boundary.
    boundary.discretise();

    // Boundary should cut diagonally through the upper and lower triangles of the element.

    // Check connectivity between nodes and boundary points.

    // First node.
    check((mesh.nodes[0].nBoundaryPoints == 2), "Incorrect number of boundary points associated with node 0!");
    check((mesh.nodes[0].boundaryPoints[0] == 0), "Boundary point 0 of node 0 is incorrect!");
    check((mesh.nodes[0].boundaryPoints[1] == 3), "Boundary point 1 of node 0 is incorrect!");

    // Second node.
    check((mesh.nodes[1].nBoundaryPoints == 2), "Incorrect number of boundary points associated with node 1!");
    check((mesh.nodes[1].boundaryPoints[0] == 0), "Boundary point 0 of node 1 is incorrect!");
    check((mesh.nodes[1].boundaryPoints[1] == 1), "Boundary point 1 of node 1 is incorrect!");

    // Third node.
    check((mesh.nodes[2].nBoundaryPoints == 2), "Incorrect number of boundary points associated with node 2!");
    check((mesh.nodes[2].boundaryPoints[0] == 2), "Boundary point 0 of node 2 is incorrect!");
    check((mesh.nodes[2].boundaryPoints[1] == 3), "Boundary point 1 of node 2 is incorrect!");

    // Fourth node.
    check((mesh.nodes[3].nBoundaryPoints == 2), "Incorrect number of boundary points associated with node 3!");
    check((mesh.nodes[3].boundaryPoints[0] == 1), "Boundary point 0 of node 3 is incorrect!");
    check((mesh.nodes[3].boundaryPoints[1] == 2), "Boundary point 1 of node 3 is incorrect!");

    // Check boundary segment connectivity.
    check((mesh.elements[0].nBoundarySegments == 2), "Incorrect number of boundary segments for node 0!");
    check((mesh.elements[0].boundarySegments[0] == 0), "Boundary segment 0 of element 0 is incorrect!");
    check((mesh.elements[0].boundarySegments[1] == 1), "Boundary segment 1 of element 0 is incorrect!");

    return 0;

error:
    return 1;
}

int testAreaFraction()
{
    // Initialise a 1x1 non-periodic mesh.
    Mesh mesh(1, 1, false);

    // Push hole into a vector container.
    std::vector<Hole> holes;
    holes.push_back(Hole(1, 1, 1));

    // Initialise the level set object.
    LevelSet levelSet(mesh, 3, holes);

    // Initialise the boundary object.
    Boundary boundary(mesh, levelSet);

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
    boundary.discretise();

    // Calculate the element area.
    boundary.computeAreaFractions();

    // Check that the element is half filled.
    check((mesh.elements[0].area == 0.5), "Element area fraction is incorrect!");

    // Sub test 2:
    // Vertical boundary cutting through the middle of the element.

    // Place nodes on left edge of element outside zero contour.
    levelSet.signedDistance[0] = -1;
    levelSet.signedDistance[2] = -1;

    // Place nodes on right edge of element inside zero contour.
    levelSet.signedDistance[1] = 1;
    levelSet.signedDistance[3] = 1;

    // Discretise the boundary.
    boundary.discretise();

    // Calculate the element area.
    boundary.computeAreaFractions();

    // Check that the element is half filled.
    check((mesh.elements[0].area == 0.5), "Element area fraction is incorrect!");

    // Sub test 3:
    // Zero contour along diagonal.

    // Place bottom left and top right nodes on zero contour.
    levelSet.signedDistance[0] = 0;
    levelSet.signedDistance[3] = 0;

    // Place bottom right node inside and top left node outside.
    levelSet.signedDistance[1] = 1;
    levelSet.signedDistance[2] = -1;

    // Discretise the boundary.
    boundary.discretise();

    // Calculate the element area.
    boundary.computeAreaFractions();

    // Check that the element is half filled.
    check((mesh.elements[0].area == 0.5), "Element area fraction is incorrect!");

    // Sub test 4:
    // Zero contour along diagonal in top right quadrant.

    // Place bottom left, bottom right, and top left nodes inside.
    levelSet.signedDistance[0] = 1;
    levelSet.signedDistance[1] = 1;
    levelSet.signedDistance[2] = 1;

    // Place top right node outside.
    levelSet.signedDistance[3] = -1;

    // Discretise the boundary.
    boundary.discretise();

    // Calculate the element area.
    boundary.computeAreaFractions();

    // Check the element area.
    // Full element minus half a quarter element, i.e. 1.0 - 1/8 = 0.875.
    check((mesh.elements[0].area == 0.875), "Element area fraction is incorrect!");

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
    boundary.discretise();

    // Calculate the element area.
    boundary.computeAreaFractions();

    // Check the element area.
    check((mesh.elements[0].area == (1.0/8.0)), "Element area fraction is incorrect!");

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
