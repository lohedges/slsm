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

#include "Optimise.h"
#include <iostream>

Optimise::Optimise(const std::vector<BoundaryPoint>& boundaryPoints_,
                   const std::vector<double>& constraintDistances_,
                   std::vector<double>& lambdas_,
                   std::vector<double>& velocities_) :
                   boundaryPoints(boundaryPoints_),
                   constraintDistances(constraintDistances_),
                   lambdas(lambdas_),
                   velocities(velocities_)
{
    nPoints = boundaryPoints.size();
    nConstraints = constraintDistances.size();
    isSideLimit.resize(nPoints);

    // Test that NlOpt wrapper function is working.

    NLOptWrapper wrapper;
    wrapper.index = 0;
    wrapper.callback = this;

    void* data = reinterpret_cast<void*>(&wrapper);

    std::vector<double> tmp1(2);
    std::vector<double> tmp2(2);

    std::cout << callbackWrapper(tmp1, tmp2, data) << '\n';
}

double Optimise::callback(const std::vector<double>& lambda, std::vector<double>& gradient, unsigned int index)
{
    // Calculate the boundary movement vector (velocities).
    computeVelocities(lambda);

    // Compute the gradients.
    computeGradients(lambda, gradient, index);

    // Return the function value.
    return computeFunction(index);
}

double callbackWrapper(const std::vector<double>& lambda, std::vector<double>& gradient, void* data)
{
    NLOptWrapper* wrapperData = reinterpret_cast<NLOptWrapper*>(data);
    return reinterpret_cast<Optimise*>(wrapperData->callback)->callback(lambda, gradient, wrapperData->index);
}

void Optimise::computeVelocities(const std::vector<double>& lambda)
{
    // Loop over all boundary points.
    for (unsigned int i=0;i<nPoints;i++)
    {
        // Reset side limit flag.
        isSideLimit[i] = false;

        // Initialise component for objective.
        velocities[i] = lambda[0]*boundaryPoints[i].sensitivities[0];

        // Add components for constraints.
        for (unsigned int j=1;j<nConstraints+1;j++)
        {
            velocities[i] += lambda[j]*boundaryPoints[i].sensitivities[j];
        }

        // Check side limits if point lies close to domain boundary.
        if (boundaryPoints[i].isDomain)
        {
            // Apply side limit.
            if (velocities[i] < boundaryPoints[i].negativeLimit)
            {
                velocities[i] = boundaryPoints[i].negativeLimit;
                isSideLimit[i] = true;
            }
        }
    }
}

double Optimise::computeFunction(unsigned int index)
{
    // This method assumes that velocities have already been calculated.

    // Intialise function.
    double func = 0;

    // Integrate function over boundary points.
    for (unsigned int i=0;i<nPoints;i++)
        func += (velocities[i] * boundaryPoints[i].sensitivities[index] * boundaryPoints[i].length);

    if (index == 0) return func;
    else return func - constraintDistances[index - 1];
}

void Optimise::computeGradients(const std::vector<double>& lambda, std::vector<double>& gradient, unsigned int index)
{
    // Whether we're at the origin, i.e. all lambda values are zero.
    bool isOrigin = false;

    // Zero the gradients.
    gradient[0] = 0;
    for (unsigned int i=1;i<nConstraints+1;i++) gradient[i] = 0;

    // If all lambda values are zero then we need a fix for the
    // analytic gradient calculation. This is because side constraints
    // for boundary points lying exactly on the domain boundary will
    // be active for lambda < 0 and inactive for lambda > 0. As such,
    // the gradient will be half as large, i.e. only positive lambda
    // contributes.

    // Sum lambda values.
    double lambdaSum = 0;
    for (unsigned int i=0;i<nConstraints+1;i++) lambdaSum += std::abs(lambda[i]);

    // If all lambdas are zero, then we're at the origin.
    if (lambdaSum < 1e-6)
    {
        isOrigin = true;

        // Find points lying on domain boundary.
        for (unsigned int i=0;i<nPoints;i++)
        {
            if (boundaryPoints[i].isDomain) isSideLimit[i] = true;
            else isSideLimit[i] = false;
        }
    }

    // Calculate the deriviative with respect to each lambda.

    // Loop over all points.
    for (unsigned int i=0;i<nPoints;i++)
    {
        // Loop over all functions (objective, then constraints).
        for (unsigned int j=0;j<nConstraints+1;j++)
        {
            if (!isSideLimit[i])
            {
                gradient[j] += (boundaryPoints[i].sensitivities[index]
                            * boundaryPoints[i].sensitivities[j]
                            * boundaryPoints[i].length);
            }
            else if (isOrigin)
            {
                gradient[j] += 0.5 * (boundaryPoints[i].sensitivities[index]
                            * boundaryPoints[i].sensitivities[j]
                            * boundaryPoints[i].length);

            }
        }
    }
}
