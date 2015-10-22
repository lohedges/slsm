/*
 * File:   level_set_tests.cpp
 * Author: lester
 */

#include "lsm.h"

int testSignedDistance()
{
    // A test for correct initialisation of the signed distance function.

    // Initialise a 10x10 non-periodic mesh.
    Mesh mesh(10, 10, false);

    // Initialise a hole.
    Hole hole;

    // Place hole in the centre of the mesh.
    hole.coord.x = 5;
    hole.coord.y = 5;
    hole.r = 3;

    // Push hole into a vector container.
    std::vector<Hole> holes;
    holes.push_back(hole);

    // Initialise the level set object.
    LevelSet levelSet(mesh, 3, holes);

    // Diagonal of reduced 3x3 box.
    double diag = sqrt(18);

    // Difference between diag and hole radius.
    double d = diag - hole.r;

    // Check signed distance against expected values.

    // Set error number.
    errno = 0;

    // Check that corner nodes lie on outer domain boundary.
    check((levelSet.signedDistance[0] == 0), "Signed distance mismatch!");
    check((levelSet.signedDistance[10] == 0), "Signed distance mismatch!");
    check((levelSet.signedDistance[110] == 0), "Signed distance mismatch!");
    check((levelSet.signedDistance[120] == 0), "Signed distance mismatch!");

    // Check points on hole surface.
    check((levelSet.signedDistance[57] == 0), "Signed distance mismatch!");
    check((levelSet.signedDistance[63] == 0), "Signed distance mismatch!");
    check((levelSet.signedDistance[27] == 0), "Signed distance mismatch!");
    check((levelSet.signedDistance[93] == 0), "Signed distance mismatch!");

    // Check points between hole surface and outer domain (inside structure).
    check((levelSet.signedDistance[56] == 1), "Signed distance mismatch!");
    check((levelSet.signedDistance[64] == 1), "Signed distance mismatch!");
    check((levelSet.signedDistance[16] == 1), "Signed distance mismatch!");
    check((levelSet.signedDistance[104] == 1), "Signed distance mismatch!");

    // Check points inside hole (outside structure).
    check((levelSet.signedDistance[58] == -1), "Signed distance mismatch!");
    check((levelSet.signedDistance[62] == -1), "Signed distance mismatch!");
    check((levelSet.signedDistance[38] == -1), "Signed distance mismatch!");
    check((levelSet.signedDistance[82] == -1), "Signed distance mismatch!");

    // Check points two nodes diagonally inside the mesh corners.
    // These should lie inside the structure by a distance equal to d.
    check((levelSet.signedDistance[24] == d), "Signed distance mismatch!");
    check((levelSet.signedDistance[30] == d), "Signed distance mismatch!");
    check((levelSet.signedDistance[90] == d), "Signed distance mismatch!");
    check((levelSet.signedDistance[96] == d), "Signed distance mismatch!");

    return 0;

error:
    return 1;
}

int all_tests()
{
    mu_suite_start();

    mu_run_test(testSignedDistance);

    return 0;
}

RUN_TESTS(all_tests);
