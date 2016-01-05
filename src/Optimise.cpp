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

double callbackWrapper(const std::vector<double>& lambda, std::vector<double>& gradient, void* data)
{
    NLoptWrapper* wrapperData = reinterpret_cast<NLoptWrapper*>(data);
    return reinterpret_cast<Optimise*>(wrapperData->callback)->callback(lambda, gradient, wrapperData->index);
}

Optimise::Optimise(const std::vector<BoundaryPoint>& boundaryPoints_,
                   const std::vector<double>& constraintDistances_,
                   std::vector<double>& lambdas_,
                   std::vector<double>& velocities_) :
                   boundaryPoints(boundaryPoints_),
                   constraintDistances(constraintDistances_),
                   lambdas(lambdas_),
                   velocities(velocities_)
{
    // Set number of points and constraints.
    nPoints = boundaryPoints.size();
    nConstraints = constraintDistances.size();

    // Resize data structures.
    isSideLimit.resize(nPoints);
    negativeLambdaLimits.resize(nConstraints + 1);
    positiveLambdaLimits.resize(nConstraints + 1);
    scaleFactors.resize(nConstraints + 1);
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

double Optimise::solve()
{
    // Compute the scale factors for the objective end constraints.
    computeScaleFactors();

    // Compute the lambda limits.
    computeLambdaLimits();

    // Create wrapper for objective.
    NLoptWrapper objectiveWrapper;
    objectiveWrapper.index = 0;
    objectiveWrapper.callback = this;

    // Create wrappers for constraints.
    NLoptWrapper constraintWrappers[nConstraints];
    for (unsigned int i=0;i<nConstraints;i++)
    {
        constraintWrappers[i].index = i + 1;
        constraintWrappers[i].callback = this;
    }

    // Instantiate NLopt optimisation object.
    nlopt::opt opt(nlopt::LD_SLSQP, 1 + nConstraints);

    // Set limits.
    opt.set_lower_bounds(negativeLambdaLimits);
    opt.set_upper_bounds(negativeLambdaLimits);

    // Set a maximum of 500 iterations.
    opt.set_maxeval(500);

    // Set convergence tolerance.
    opt.set_xtol_rel(1e-6);

    // Specify that we want to minimise the objective function.
    opt.set_min_objective(callbackWrapper, &objectiveWrapper);

    // Add the inequality constraints.
    for (unsigned int i=0;i<nConstraints;i++)
        opt.add_inequality_constraint(callbackWrapper, &constraintWrappers[i], 1e-8);

    // The optimum value of the objective function.
    double optObjective;

    // Perform the optimisation.
    returnCode = opt.optimize(lambdas, optObjective);

    return optObjective;
}

void Optimise::computeScaleFactors()
{

}

void Optimise::computeLambdaLimits()
{

}

void Optimise::computeVelocities(const std::vector<double>& lambda)
{
    // Loop over all boundary points.
    for (unsigned int i=0;i<nPoints;i++)
    {
        // Reset side limit flag.
        isSideLimit[i] = false;

        // Initialise component for objective.
        velocities[i] = scaleFactors[0] * lambda[0] * boundaryPoints[i].sensitivities[0];

        // Add components for constraints.
        for (unsigned int j=1;j<nConstraints+1;j++)
        {
            velocities[i] += scaleFactors[j] * lambda[j] * boundaryPoints[i].sensitivities[j];
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
        func += (scaleFactors[index] * velocities[i] * boundaryPoints[i].sensitivities[index] * boundaryPoints[i].length);

    if (index == 0) return func;
    else return func - constraintDistances[index - 1];
}

void Optimise::computeGradients(const std::vector<double>& lambda, std::vector<double>& gradient, unsigned int index)
{
    // Whether we're at the origin, i.e. all lambda values are zero.
    bool isOrigin = false;

    // Scale factor.
    double scaleFactor = scaleFactors[index] * scaleFactors[index];

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

    // Calculate the derivative with respect to each lambda.

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
                            * boundaryPoints[i].length
                            * scaleFactor);
            }
            else if (isOrigin)
            {
                gradient[j] += (boundaryPoints[i].sensitivities[index]
                            * boundaryPoints[i].sensitivities[j]
                            * boundaryPoints[i].length
                            * 0.5 * scaleFactor);
            }
        }
    }
}

void Optimise::queryReturnCode()
{
    // Success.
    if (returnCode == 1) std::cout << "[INFO] Success: Generic success return value.\n";
    else if (returnCode == 2) std::cout << "[INFO] Success: stopval was reached.\n";
    else if (returnCode == 3) std::cout << "[INFO] Success: ftol_rel or ftol_abs was reached.\n";
    else if (returnCode == 4) std::cout << "[INFO] Success: xtol_rel or xtol_abs was reached.\n";
    else if (returnCode == 5) std::cout << "[INFO] Success: maxeval was reached.\n";
    else if (returnCode == 6) std::cout << "[INFO] Success: maxtime was reached.\n";

    // Error codes.
    else if (returnCode == -1) std::cout << "[INFO] Failed: Generic failure code.\n";
    else if (returnCode == -2) std::cout << "[INFO] Failed: Invalid arguments.\n";
    else if (returnCode == -3) std::cout << "[INFO] Failed: Ran out of memory.\n";
    else if (returnCode == -4) std::cout << "[INFO] Failed: Halted because roundoff errors limited progress.\n";
    else if (returnCode == -5) std::cout << "[INFO] Failed: Halted because of a forced termination.\n";

    // No optimisation.
    else std::cout << "[INFO] No optimisation has been performed.\n";
}
