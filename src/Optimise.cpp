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
    nConstraints = lambdas.size() - 1;

    // Check for empty constraint distances vector.
    if (constraintDistances.empty())
    {
        errno = 0;
        log_warn("Empty constraint distances vector. Assuming unconstrained optimisation.");
        nConstraints = 0;
    }

    // Resize data structures.
    isSideLimit.resize(nPoints);
    negativeLambdaLimits.resize(nConstraints + 1);
    positiveLambdaLimits.resize(nConstraints + 1);
    scaleFactors.resize(nConstraints + 1);

    // Copy constraint distances.
    constraintDistancesScaled = constraintDistances;
}

double Optimise::callback(const std::vector<double>& lambda, std::vector<double>& gradient, unsigned int index)
{
    // Calculate the boundary movement vector (velocities).
    computeVelocities(lambda);

    // Compute the gradients.
    if (!gradient.empty())
        computeGradients(lambda, gradient, index);

    // Return the function value.
    return computeFunction(index);
}

double Optimise::solve()
{
    // Compute the scale factors for the objective end constraints.
    computeScaleFactors();

    // Compute scaled constraint change distances.
    computeConstraintDistances();

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
    opt.set_upper_bounds(positiveLambdaLimits);

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

    // Unscale the objective function change.
    optObjective /= scaleFactors[0];

    // Compute the optimum velocity vector.
    computeVelocities(lambdas);

    // Rescale the velocities and lambda values (if necessary).
    rescaleVelocities();

    // Unscale the lambda values.
    for (unsigned int i=0;i<nConstraints+1;i++)
        lambdas[i] *= scaleFactors[i];

    return optObjective;
}

void Optimise::computeScaleFactors()
{
    /* In order for the optimiser to work effectively it is important
       that small changes in the lambdas result in small changes in the
       functions for the objective and constraints and their respective
       gradients. Since we are solving a multi-dimensional optimisation
       problem it is important that all variables are on the same scale.
       This enables us to use a universal convergence tolerance.

       Our scaling protocol is described below. Note that we choose to
       store scale factors for each function, rather than scaling (then
       rescaling) the input data. Sensitivites are scaled down (reduced)
       and the lambda values are scale up (enlarged).

         1) For each function, scale by the largest absolute sensitivity,
            i.e. the maximum magnitude is one.

         2) Scale by the gradient at the origin (all lambdas are zero).
            This ensures that the gradient is close to one.

        There is no need to independently scale the boundary integral
        coefficients since the boundary lengths are independent of the
        function, i.e.

            c^f_i = s^f_i * l_i

        The l_i variables are the same for the objective and constraints
        so, by definition, the c^f's are on the same scale if we simply
        normalise by max(abs(s^f_i)).
     */

    // Loop over all functions: objective first, then constraints.
    for (unsigned int i=0;i<nConstraints+1;i++)
    {
        // Initialise maximum sensitivity.
        double maxSens = 0;

        // Loop over all boundary points.
        for (unsigned int j=0;j<nPoints;j++)
        {
            // Test whether sensitivity magnitude is current maximum.
            double sens = std::abs(boundaryPoints[j].sensitivities[i]);
            if (sens > maxSens) maxSens = sens;
        }

        // Store scale factor.
        scaleFactors[i] = (1.0 / maxSens);
    }

    // Create lambda vector (all zeros, i.e. at the origin).
    std::vector<double> lambda(nConstraints + 1);
    std::fill(lambda.begin(), lambda.end(), 0.0);

    // Initialise gradient vector.
    std::vector<double> gradient(nConstraints + 1);

    // Loop over all functions: objective first, then constraints.
    for (unsigned int i=0;i<nConstraints+1;i++)
    {
        // Calculate the gradient.
        computeGradients(lambda, gradient, i);

        // Scale by diagonal gradient entry (absolute value).
        scaleFactors[i] *= (1.0 / std::abs(gradient[i]));
    }
}

void Optimise::computeConstraintDistances()
{
    /* If we are far from satisfying the constraint then we need
       to scale the constraint distance so that it can be "satisfied"
       by simply moving in the correct direction, i.e. moving towards
       satisying the constraint.

       I would have thought that the optimiser should be able to minimise
       the constraint violation (this is suggested in the SLP paper) but
       this doesn't seem to be the case, i.e. the optimiser would return
       the lambdas that minimise the objective and the constraint violation
       (if the constraint can't be satisfied).
     */

    // Loop over all constraints.
    for (unsigned int i=0;i<nConstraints;i++)
    {
        // Min and max changes.
        double min = 0;
        double max = 0;

        // Integrate over boundary points.
        for (unsigned int j=0;j<nPoints;j++)
        {
            if (boundaryPoints[j].sensitivities[i+1] > 0)
            {
                min += boundaryPoints[j].sensitivities[i+1]
                     * boundaryPoints[j].length
                     * boundaryPoints[j].negativeLimit;

                max += boundaryPoints[j].sensitivities[i+1]
                     * boundaryPoints[j].length
                     * boundaryPoints[j].positiveLimit;
            }
            else
            {
                min += boundaryPoints[j].sensitivities[i+1]
                     * boundaryPoints[j].length
                     * boundaryPoints[j].positiveLimit;

                max += boundaryPoints[j].sensitivities[i+1]
                     * boundaryPoints[j].length
                     * boundaryPoints[j].negativeLimit;
            }
        }

        // Scale (20% is arbitrary, but seems to work well).
        min *= 0.2;
        max *= 0.2;

        // We want to reduce the constraint function.
        if (constraintDistances[i] < 0)
        {
            if (constraintDistances[i] < min)
                constraintDistancesScaled[i] = min;
        }
        // Raise constraint function.
        else
        {
            if (constraintDistances[i] > max)
                constraintDistancesScaled[i] = max;
        }
    }
}

void Optimise::computeLambdaLimits()
{
    // Dummy values for now.
    for (unsigned int i=0;i<nConstraints+1;i++)
    {
        negativeLambdaLimits[i] = -200;
        positiveLambdaLimits[i] = 200;
    }
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
            velocities[i] += scaleFactors[j] * lambda[j] * boundaryPoints[i].sensitivities[j];

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
    else return (func - (scaleFactors[index] * constraintDistancesScaled[index - 1]));
}

void Optimise::computeGradients(const std::vector<double>& lambda, std::vector<double>& gradient, unsigned int index)
{
    // Whether we're at the origin, i.e. all lambda values are zero.
    bool isOrigin = false;

    // Zero the gradients.
    gradient[0] = 0;
    for (unsigned int i=1;i<nConstraints+1;i++) gradient[i] = 0;

    /* If all lambda values are zero then we need a fix for the
       analytic gradient calculation. This is because side constraints
       for boundary points lying exactly on the domain boundary will
       be active for lambda < 0 and inactive for lambda > 0. As such,
       the gradient will be half as large, i.e. only positive lambda
       contributes.
     */

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
            // Scale factor.
            double scaleFactor = scaleFactors[index] * scaleFactors[j];

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

void Optimise::rescaleVelocities()
{
    // Check for CFL violation and rescale the velocities
    // and lambda values if necessary.

    // Maximum velocity magnitude.
    double maxVel = 0;

    // Velocity to rescale to (site dependent CFL limit).
    double vRescale;

    // Loop over all boundary points.
    for (unsigned int i=0;i<nPoints;i++)
    {
        // Whether CFL is violated.
        bool isCFL = false;

        if (velocities[i] < boundaryPoints[i].negativeLimit)
        {
            vRescale = -boundaryPoints[i].negativeLimit;
            isCFL = true;
        }
        else if (velocities[i] > boundaryPoints[i].positiveLimit)
        {
            vRescale = boundaryPoints[i].positiveLimit;
            isCFL = true;
        }

        if (isCFL)
        {
            double vel = std::abs(velocities[i]);

            // Check if current maximum is exceeded.
            if (vel > maxVel) maxVel = vel;
        }
    }

    // CFL condition is violated, rescale velocities.
    if (maxVel)
    {
        // Scaling factor.
        double scale = vRescale / maxVel;

        // Scale lambda values.
        for (unsigned int i=0;i<nConstraints+1;i++)
            lambdas[i] *= scale;

        // Recompute the velocity vector.
        computeVelocities(lambdas);
    }
}

void Optimise::queryReturnCode()
{
    /* N.B.
        Despite providing an extensive list of return codes NLopt does not
        report when the optmisation exits because a lower or upper lambda
        limit is hit. In this case the return code is 4. However, it's easy
        to test for this situation by comparing the "optimum" lambda values
        to the limits.
     */

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
