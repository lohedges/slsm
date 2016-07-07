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

int testUpwindFiniteDifference()
{
    // A test for a correct upwind finite difference
    // approximation for negative phi. This was a bug
    // in scikit-fmm.

    // Initialise a 4x4 non-periodic mesh.
    slsm::Mesh mesh(4, 4, false);

    // Create a single hole (just to pass to constructor).
    std::vector<slsm::Hole> hole;

    // Initialise the level set object.
    slsm::LevelSet levelSet(mesh, hole, 0.5, 3);

    // Place all nodes outside the structure.
    for (unsigned int i=0;i<mesh.nNodes;i++)
        levelSet.signedDistance[i] = -1;

    // Place nodes at bottom left and top right of mesh inside the structure.
    levelSet.signedDistance[0] = 1;
    levelSet.signedDistance[mesh.nNodes-1] = 1;

    // Fill array with expected values.
    double expected[25] =
        {0.35355339, -0.50000000, -1.45118446, -2.43491262, -3.23422652,
        -0.50000000, -1.20710678, -2.00007282, -2.73579936, -2.43491262,
        -1.45118446, -2.00007282, -2.65444013, -2.00007282, -1.45118446,
        -2.43491262, -2.73579936, -2.00007282, -1.20710678, -0.50000000,
        -3.23422652, -2.43491262, -1.45118446, -0.50000000,  0.35355339};

    // Reinitialise the signed distance function.
    levelSet.reinitialise();

    // Set error number.
    errno = 0;

    // Check signed distance against expected values.
    for (unsigned int i=0;i<mesh.nNodes;i++)
        lsm_check((std::abs(levelSet.signedDistance[i] - expected[i]) < 1e-6), "Signed distance mismatch!");

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
