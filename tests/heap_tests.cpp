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

#include "lsm.h"

/* N.B.
    The Heap class can perform self testing when passed isTest = true
    as a constructor argument. These tests check basic functionality
    of the heap beyond what is performed within the class itself.
 */

#include "lsm.h"

int testPush()
{
    // Initialise heap.
    Heap heap(10, true);

    // Set error number.
    errno = 0;

    // Add entries to the heap.
    for (int i=0;i<10;i++)
    {
        double value = -i;

        // Push onto the heap.
        heap.push(i, value);

        // Check that heap size is correct.
        check(heap.size() == unsigned(i + 1), "Heap push: incorrect heap size!");

        // Check that new value has sifted to the top of the heap.
        check(heap.peek() == value, "Heap push: incorrect value at top of heap!");
    }

    return 0;

error:
    return 1;
}

int testPop()
{
    // Initialise random number generator.
    MersenneTwister rng;

    // Initialise vector of doubles.
    std::vector<double> vec(10);

    // Fill vector with random numbers.
    for (unsigned int i=0;i<vec.size();i++)
        vec[i] = rng();

    // Copy vector.
    std::vector<double> sorted = vec;

    // Sort vector.
    std::sort(sorted.begin(), sorted.end());

    // Initialise heap.
    Heap heap(vec.size(), true);

    // Initialise back pointer array.
    std::vector<unsigned int> heapPtr(vec.size());

    // Push values onto the heap.
    for (unsigned int i=0;i<vec.size();i++)
        heapPtr[i] = heap.push(i, vec[i]);

    // Set error number.
    errno = 0;

    // Now pop off values.
    for (unsigned int i=0;i<vec.size();i++)
    {
        unsigned int addr;
        double value;

        heap.pop(addr, value);

        // Make sure that values are in ascending order.
        check(value == sorted[i], "Heap pop: incorrect value!");
    }

    return 0;

error:
    return 1;
}

int testSet()
{
    // Initialise random number generator.
    MersenneTwister rng;

    // Initialise vector of doubles.
    std::vector<double> vec(10);

    // Fill vector with random numbers.
    for (unsigned int i=0;i<vec.size();i++)
        vec[i] = rng();

    // Copy vector.
    std::vector<double> sorted = vec;

    // Sort vector.
    std::sort(sorted.begin(), sorted.end());

    // Initialise heap.
    Heap heap(vec.size(), true);

    // Initialise back pointer array.
    std::vector<unsigned int> heapPtr(vec.size());

    // Push values onto the heap.
    for (unsigned int i=0;i<vec.size();i++)
        heapPtr[i] = heap.push(i, vec[i]);

    // Set 5th value to large negative number.
    heap.set(heapPtr[4], -100);

    // Set error number.
    errno = 0;

    // Check that top entry is correct.
    check(heap.peek() == -100, "Heap set: incorrect value!");

    return 0;

error:
    return 1;
}

int all_tests()
{
    mu_suite_start();

    mu_run_test(testPush);
    mu_run_test(testPop);
    mu_run_test(testSet);

    return 0;
}

RUN_TESTS(all_tests);
