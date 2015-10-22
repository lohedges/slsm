/*
 * File:   fast_marching_tests.cpp
 * Author: lester
 */

#include "lsm.h"

int testUpwindFiniteDifference()
{
    // A test for a correct upwind finite difference
    // approximation for negative phi. This was a bug
    // in scikit-fmm.

    // Initialise a 4x4 non-periodic mesh.
    Mesh mesh(4, 4, false);

    // Create a single hole (just to pass to constructor).
    std::vector<Hole> hole;

    // Initialise the level set object.
    LevelSet levelSet(mesh, 3, hole);

    // Place all nodes outside the structure.
    for (unsigned int i=0;i<levelSet.nNodes;i++)
        levelSet.signedDistance[i] = -1;

    // Place nodes at bottom left and top right of mesh inside the structure.
    levelSet.signedDistance[0] = 1;
    levelSet.signedDistance[levelSet.nNodes-1] = 1;

    // Fill array with expected values.
    double expected[25] =
        {0.35355339, -0.50000000, -1.45118446, -2.43491262, -3.23422652,
        -0.50000000, -1.20710678, -2.00007282, -2.73579936, -2.43491262,
        -1.45118446, -2.00007282, -2.65444013, -2.00007282, -1.45118446,
        -2.43491262, -2.73579936, -2.00007282, -1.20710678, -0.50000000,
        -3.23422652, -2.43491262, -1.45118446, -0.50000000,  0.35355339};

    // Re-initialise the signed distance function.
    levelSet.reinitialise();

    // Set error number.
    errno = 0;

    // Check signed distance against expected values.
    for (unsigned int i=0;i<levelSet.nNodes;i++)
        check((std::abs(levelSet.signedDistance[i] - expected[i]) < 1e-6), "Signed distance mismatch!");

    return 0;

error:
    return 1;
}

int all_tests()
{
    mu_suite_start();

    mu_run_test(testUpwindFiniteDifference);

    return 0;
}

RUN_TESTS(all_tests);
