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

#ifndef _OPTIMISE_H
#define _OPTIMISE_H

#include <iostream>
#include <nlopt.hpp>

#include "Boundary.h"

// ASSOCIATED DATA TYPES

//! Wrapper structure for interfacing with NLopt.
struct NLoptWrapper
{
    unsigned int index;     //!< Function index (0 = objective, 1, 2, ... = constraints).
    void* callback;         //!< Pointer to callback function wrapper.
};

//! NLopt callback function wrapper.
/*! \param lambda
        The current lambda value for each function.

    \param gradient
        The gradient of the change in objective or constraint with respect to each lambda.

    \param data
        Void pointer to NLoptWrapper data.
 */
double callbackWrapper(const std::vector<double>& lambda, std::vector<double>& gradient, void* data);

/*!\brief A class to solve for the optimum boundary movement vector (velocity).

    We make use of the NLopt SLSQP solver:
        http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms#SLSQP
 */
class Optimise
{
public:
    //! Constructor.
    /*! \param boundaryPoints_
            A reference to a vector of boundary points.

        \param constraintDistances_
            Distance from violation for each constraint (negative values
            indicate that we want to reduce the constraint function).

        \param lambdas_
            The optimum lambda values. This array is modified.

        \param velocities_
            The optimum boundary movement vector (velocity). This array is modified.
     */
    Optimise(const std::vector<BoundaryPoint>&, const std::vector<double>&,
        std::vector<double>&, std::vector<double>&);

    //! Execute the NLopt SLSQP solver.
    /*! \return
            The optimum value of the objective function.
     */
    double solve();

    //! NLopt callback function.
    /*! \param lambda
            The current lambda value for each function.

        \param gradient
            The gradient of the change in objective or constraint with respect to each lambda.

        \param index
            Function index, 0 = objective, 1, 2, 3, ... = constraints.
     */
    double callback(const std::vector<double>&, std::vector<double>&, unsigned int);

    //! Query the NLopt return code.
    void queryReturnCode();

private:
    /// The number of boundary points.
    unsigned int nPoints;

    /// The number of constraints.
    unsigned int nConstraints;

    /// A reference to a vector of boundary points.
    const std::vector<BoundaryPoint>& boundaryPoints;

    /// A reference to a vector of constraint violation distances.
    const std::vector<double>& constraintDistances;

    /// A reference to a vector of optimum lambda values (to be found by solver).
    std::vector<double>& lambdas;

    /// A reference to a vector of optimum velocity values (to be found by solver).
    std::vector<double>& velocities;

    /// Whether the velocity side limit is active for each boundary point.
    std::vector<bool> isSideLimit;

    /// Negative lambda limits.
    std::vector<double> negativeLambdaLimits;

    /// Positive lambda limits.
    std::vector<double> positiveLambdaLimits;

    /// Scale factor for each function.
    std::vector<double> scaleFactors;

    /// Scaled constraint distances.
    std::vector<double> constraintDistancesScaled;

    /// Optimiser return code.
    nlopt::result returnCode;

    //! Compute scale factors.
    void computeScaleFactors();

    //! Compute constraint distances.
    void computeConstraintDistances();

    //! Compute lambda limits.
    void computeLambdaLimits();

    //! Compute the boundary movement vector (velocities).
    /*! \param lambda
            A vector of lambda values (objective, then constraints).
     */
    void computeVelocities(const std::vector<double>&);

    //! Compute the change in the objective or constraint functions.
    /*! \param index
            Function index, 0 = objective, 1, 2, 3, ... = constraints.
     */
    double computeFunction(unsigned int);

    //! Compute the gradient of the objective or constraint functions.
    /*! \param lambda
            A vector of lambda values (objective, then constraints).

        \param gradient
            The gradient of the change in objective or constraint with respect to each lambda.

        \param index
            Function index, 0 = objective, 1, 2, 3, ... = constraints.
     */
    void computeGradients(const std::vector<double>&, std::vector<double>&, unsigned int index);

    //! Rescale velocities and lambda values if the CFL condition is violated.
    void rescaleVelocities();
};

#endif  /* _OPTIMISE_H */
