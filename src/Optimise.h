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

/*! \file Optimise.h
    \brief A class for finding the solution for the optimum velocity vector.
 */

#include <iostream>
#include <limits>
#include <numeric>
#include <nlopt.hpp>

#include "Boundary.h"
#include "Debug.h"

namespace lsm
{
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

    /*!\brief A class to solve for the optimum velocity vector.

        Find the velocity vector that minimises (or maximises) the change in the
        objective function subject to an arbitrary number of inequality constraints.

        Support for equality constraints will be added in a future version.

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
                Distance from each constraint (negative values indicate that
                the constraint is satisfied).

            \param lambdas_
                The optimum lambda values. This array is modified.

            \param timeStep_
                The effective time step. The optimum boundary displacement vector is the
                velocity vector multiplied by the time step. In many level set problems
                the time step is assumed to be one, i.e. the displacement and velocity
                vectors are equivalent. The effective time step is taken as the absolute
                value of the lambda value for the objective, i.e. abs(lambdas[0]).

            \param maxDisplacement_
                The maximum displacement (default is 0.5).

            \param isMax_
                Whether to maximise the objective function (default is to minimise).
         */
        Optimise(std::vector<BoundaryPoint>&, const std::vector<double>&,
            std::vector<double>&, double&, double maxDisplacement_ = 0.5, bool isMax_ = false);

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

        /// The number of active constraints.
        unsigned int nConstraints;

        /// The number of initial constraints.
        unsigned int nConstraintsInitial;

        /// A reference to a vector of boundary points.
        std::vector<BoundaryPoint>& boundaryPoints;

        /// A reference to a vector of constraint violation distances.
        const std::vector<double>& constraintDistances;

        /// A reference to a vector of optimum lambda values (to be found by solver).
        std::vector<double>& lambdas;

        /// The effective time step.
        double& timeStep;

        /// The maximum displacement.
        double maxDisplacement;

        /// Whether to maximise the objective function.
        bool isMax;

        /// Whether constraints are active.
        std::vector<bool> isActive;

        /// A map between indices for active constraints.
        std::vector<unsigned int> indexMap;

        /// The boundary point displacement vector.
        std::vector<double> displacements;

        /// Whether the displacement side limit is active for each boundary point.
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

        //! Compute the boundary movement vector.
        /*! \param lambda
                A vector of lambda values (objective, then constraints).
         */
        void computeDisplacements(const std::vector<double>&);

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

        //! Rescale displacements and lambda values if the CFL condition is violated.
        /*! \return
                The scale factor.
         */
        double rescaleDisplacements();
    };
}

#endif  /* _OPTIMISE_H */
