/*
 * File:   lsm.cpp
 * Author: lester
 */

#include <iostream>
#include <fstream>

#include "lsm.h"

int main(int argc, char** argv)
{
    // Print git commit info, if present.
#ifdef COMMIT
    std::cout << "Git commit: " <<  COMMIT << "\n";
#endif

    // Initialise a mesh.
    Mesh mesh(160, 80, false);

    // Create a hole.
    /*std::vector<Hole> holes;
    holes.push_back(Hole(5, 5, 5));

    // Initialise the level set function (from hole vector).
    LevelSet levelSet(mesh, 3, holes);*/

    // Initialise the level set function (default Swiss cheese).
    LevelSet levelSet(mesh, 3);

    // Initialise io object.
    InputOutput io;

    // Read Peter's level set.
    std::ifstream infile("input.txt");
    unsigned int i = 0;
    while (infile >> levelSet.signedDistance[i])
        i++;

    // Re-initialise the level set to a signed distance function.
    levelSet.reinitialise();

    // Initalise the boundary object.
    Boundary boundary(mesh, levelSet);

    // Discretise the boundary.
    boundary.discretise();

    for (unsigned int i=0;i<boundary.nPoints;i++)
    {
        std::cout << boundary.points[i].x
            << ' ' << boundary.points[i].y << '\n';
    }

    // Save Paraview LSF file.
    io.saveLevelSetVTK(1, mesh, levelSet);

    return (EXIT_SUCCESS);
}
