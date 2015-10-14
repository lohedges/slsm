/*
 * File:   mesh_tests.cpp
 * Author: lester
 */

#include "lsm.h"

int testMeshSize()
{
    // Initialise a 3x3 periodic mesh.
    Mesh mesh(3, 3, true);

    check(mesh.width == 3, "Mesh width is incorrect!");
    check(mesh.height == 3, "Mesh height is incorrect!");
    check(mesh.nElements == 9, "Number of elements is incorrect!");
    check(mesh.nNodes == 16, "Number of nodes is incorrect!");

    return 0;

error:
    return 1;
}

int testNodeCoordinates()
{
    // Initialise a 3x3 periodic mesh.
    Mesh mesh(3, 3, true);

    // Check coordinates of 1st node (bottom left).
    check(mesh.nodes[0].coord.x == 0, "x coordinate of node 0 is incorrect!");
    check(mesh.nodes[0].coord.y == 0, "y coordinate of node 0 is incorrect!");

    // Check coordinates of 5th node (bulk).
    check(mesh.nodes[5].coord.x == 1, "x coordinate of node 5 is incorrect!");
    check(mesh.nodes[5].coord.y == 1, "y coordinate of node 5 is incorrect!");

    // Check coordinates of 15th node (top right).
    check(mesh.nodes[15].coord.x == 3, "x coordinate of node 5 is incorrect!");
    check(mesh.nodes[15].coord.y == 3, "y coordinate of node 5 is incorrect!");

    return 0;

error:
    return 1;
}

int testNodeConnectivity()
{
    // Initialise a 3x3 periodic mesh.
    Mesh mesh(3, 3, true);

    // Initialise a 3x3 non-periodic mesh.
    Mesh npMesh(3, 3, false);

    // Check nearest neighbours of 1st node (bottom left).
    check(mesh.nodes[0].neighbours[0] == 3, "Periodic mesh: Neighbour 0 of node 0 is incorrect!");
    check(mesh.nodes[0].neighbours[1] == 1, "Periodic mesh: Neighbour 1 of node 0 is incorrect!");
    check(mesh.nodes[0].neighbours[2] == 12, "Periodic mesh: Neighbour 2 of node 0 is incorrect!");
    check(mesh.nodes[0].neighbours[3] == 4, "Periodic mesh: Neighbour 3 of node 0 is incorrect!");

    // Check nearest neighbours of 5th node (bulk).
    check(mesh.nodes[5].neighbours[0] == 4, "Periodic mesh: Neighbour 0 of node 5 is incorrect!");
    check(mesh.nodes[5].neighbours[1] == 6, "Periodic mesh: Neighbour 1 of node 5 is incorrect!");
    check(mesh.nodes[5].neighbours[2] == 1, "Periodic mesh: Neighbour 2 of node 5 is incorrect!");
    check(mesh.nodes[5].neighbours[3] == 9, "Periodic mesh: Neighbour 3 of node 5 is incorrect!");

    // Check nearest neighbours of 15th node (top right).
    check(mesh.nodes[15].neighbours[0] == 14, "Periodic mesh: Neighbour 0 of node 15 is incorrect!");
    check(mesh.nodes[15].neighbours[1] == 12, "Periodic mesh: Neighbour 1 of node 15 is incorrect!");
    check(mesh.nodes[15].neighbours[2] == 11, "Periodic mesh: Neighbour 2 of node 15 is incorrect!");
    check(mesh.nodes[15].neighbours[3] == 3, "Periodic mesh: Neighbour 3 of node 15 is incorrect!");

    // Check nearest neighbours of 1st node (bottom left).
    check(npMesh.nodes[0].neighbours[0] == -1, "Non-periodic mesh: Neighbour 0 of node 0 is incorrect!");
    check(npMesh.nodes[0].neighbours[1] == 1, "Non-periodic mesh: Neighbour 1 of node 0 is incorrect!");
    check(npMesh.nodes[0].neighbours[2] == -1, "Non-periodic mesh: Neighbour 2 of node 0 is incorrect!");
    check(npMesh.nodes[0].neighbours[3] == 4, "Non-periodic mesh: Neighbour 3 of node 0 is incorrect!");

    // Check nearest neighbours of 5th node (bulk).
    check(npMesh.nodes[5].neighbours[0] == 4, "Non-periodic mesh: Neighbour 0 of node 5 is incorrect!");
    check(npMesh.nodes[5].neighbours[1] == 6, "Non-periodic mesh: Neighbour 1 of node 5 is incorrect!");
    check(npMesh.nodes[5].neighbours[2] == 1, "Non-periodic mesh: Neighbour 2 of node 5 is incorrect!");
    check(npMesh.nodes[5].neighbours[3] == 9, "Non-periodic mesh: Neighbour 3 of node 5 is incorrect!");

    // Check nearest neighbours of 15th node (top right).
    check(npMesh.nodes[15].neighbours[0] == 14, "Non-periodic mesh: Neighbour 0 of node 15 is incorrect!");
    check(npMesh.nodes[15].neighbours[1] == -1, "Non-periodic mesh: Neighbour 1 of node 15 is incorrect!");
    check(npMesh.nodes[15].neighbours[2] == 11, "Non-periodic mesh: Neighbour 2 of node 15 is incorrect!");
    check(npMesh.nodes[15].neighbours[3] == -1, "Non-periodic mesh: Neighbour 3 of node 15 is incorrect!");

    return 0;

error:
    return 1;
}

int testElementNodeConnectivity()
{
    // Initialise a 3x3 periodic mesh.
    Mesh mesh(3, 3, true);

    // Check the four nodes of the 1st element.
    check(mesh.elements[0].nodes[0] == 0, "Node 0 of element 0 is incorrect!");
    check(mesh.elements[0].nodes[1] == 1, "Node 1 of element 0 is incorrect!");
    check(mesh.elements[0].nodes[2] == 4, "Node 2 of element 0 is incorrect!");
    check(mesh.elements[0].nodes[3] == 5, "Node 3 of element 0 is incorrect!");

    return 0;

error:
    return 1;
}

int testNodeElementConnectivity()
{
    // Initialise a 3x3 periodic mesh.
    Mesh mesh(3, 3, true);

    // Check the elements connected to the 1st node (one element connected).
    check(mesh.nodes[0].nElements == 1, "Number of elements connected to node 0 is incorrect!");
    check(mesh.nodes[0].elements[0] == 0, "Index of element 0 connected to node 0 is incorrect!");

    // Check the elements connected to the 7th node (two elements connected).
    check(mesh.nodes[7].nElements == 2, "Number of elements connected to node 0 is incorrect!");
    check(mesh.nodes[7].elements[0] == 2, "Index of element 0 connected to node 7 is incorrect!");
    check(mesh.nodes[7].elements[1] == 5, "Index of element 0 connected to node 7 is incorrect!");

    // Check the elements connected to the 10th node (four elements connected).
    check(mesh.nodes[10].nElements == 4, "Number of elements connected to node 0 is incorrect!");
    check(mesh.nodes[10].elements[0] == 4, "Index of element 0 connected to node 7 is incorrect!");
    check(mesh.nodes[10].elements[1] == 5, "Index of element 1 connected to node 7 is incorrect!");
    check(mesh.nodes[10].elements[2] == 7, "Index of element 2 connected to node 7 is incorrect!");
    check(mesh.nodes[10].elements[3] == 8, "Index of element 3 connected to node 7 is incorrect!");

    return 0;

error:
    return 1;
}

int all_tests()
{
    mu_suite_start();

    mu_run_test(testMeshSize);
    mu_run_test(testNodeCoordinates);
    mu_run_test(testNodeConnectivity);
    mu_run_test(testElementNodeConnectivity);
    mu_run_test(testNodeElementConnectivity);

    return 0;
}

RUN_TESTS(all_tests);
