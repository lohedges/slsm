/*
 * File:   demo.cpp
 * Author: lester
 */

#include <iostream>

#include "lsm.h"

int main(int argc, char** argv)
{
    // Print git commit info, if present.
#ifdef COMMIT
    std::cout << "Git commit: " <<  COMMIT << "\n";
#endif

    // Initialise a mesh.
    Mesh mesh(10, 10, false);

    // Create a hole.
    std::vector<Hole> holes;
    holes.push_back(Hole(5, 5, 5));

    // Initialise the level set function.
    LevelSet levelSet(mesh, 3, holes);

    // Initialise a heap.
    /*Heap heap(mesh.nNodes, true);

    // Create a random vector.
    std::vector<double> vec(10);

    MersenneTwister rng;

    for (unsigned int i=0;i<vec.size();i++)
    {
        vec[i] = rng();
        std::cout << i << ' ' << vec[i] << '\n';
    }
    std::cout << '\n';

    // Push values onto heap.
    for (unsigned int i=0;i<vec.size();i++)
    {
        heap.push(i, std::abs(vec[i]));
    }

    std::cout << heap.size() << ' ' << heap.peek() << '\n';*/

    FastMarchingMethod fmm(mesh, true);

    std::vector<double> lsf(10);
    std::vector<double> vel(10);

    fmm.march(levelSet.signedDistance);

    for (unsigned int i=0;i<levelSet.nNodes;i++)
    {
        std::cout << mesh.nodes[i].coord.x << ' '
            << mesh.nodes[i].coord.y << ' '
            << levelSet.signedDistance[i] << '\n';
    }

    return (EXIT_SUCCESS);
}
