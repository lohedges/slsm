/*
  Copyright (c) 2015-2016 Lester Hedges <lester.hedges+slsm@gmail.com>

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

#include <nlopt.hpp>

namespace slsm
{
    // FORWARD DECLARATIONS

    struct BoundaryPoint;

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
        objective function subject to an arbitrary number of equality and inequality
        constraints.

        Optimisation makes use of the NLopt library:
            http://ab-initio.mit.edu/wiki/index.php/NLopt

        Note that the optimisation problem is only non-linear due to the presence
        of side constraints on the movement of boundary points lying on, or close
        to, the domain boundary. Points are not allowed to move outside the domain.

        By default, the optmiser makes use of the SLSQP algorithm, which provides
        support for equality and inequality constraints.

            http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms#SLSQP

        Support for alternative algorithms is provided through an optional constructor
        argument. Note that the chosen algorithm must support the type of constraints
        that are imposed in the optimisation problem. Check the NLopt algorithms page
        for details:

            http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms
     */
    class Optimise
    {
    public:
        //! Constructor.
        /*! \param boundaryPoints_
                A reference to a vector of boundary points.

            \param constraintDistances_
                Distance from each constraint (negative values indicate that the
                constraint is satisfied).

            \param lambdas_
                The optimum lambda values. This array is modified.

            \param timeStep_
                The effective time step. The optimum boundary displacement vector is the
                velocity vector multiplied by the time step. In many level set problems
                the time step is assumed to be one, i.e. the displacement and velocity
                vectors are equivalent. The effective time step is taken as the absolute
                value of the lambda for the objective, i.e. abs(lambdas[0]). (Note that
                for minimisation problems lambdas[0] is negative by construction so
                abs(lambdas[0]) is the same as -lambdas[0].)

            \param maxDisplacement_
                (Optional) The maximum displacement (default = 0.5).

            \param isMax_
                (Optional) Whether to maximise the objective function (default = minimise).

            \param isEquality_
                (Optional) Whether each constraint is an equality (default = inequality).

            \param algorithm_
                (Optional) The NLopt algorithm (default = LD_SLSQP).
         */
        Optimise(std::vector<BoundaryPoint>&, std::vector<double>, std::vector<double>&,
            double&, double maxDisplacement_ = 0.5, bool isMax_ = false,
            const std::vector<bool>& isEquality_ = {}, nlopt::algorithm algorithm_ = nlopt::LD_SLSQP);

        //! Execute the NLopt solver.
        /*! \return
                The optimum change in the objective function.
         */
        double solve();

        //! NLopt callback function.
        /*! \param lambda
                The current lambda value for each function.

            \param gradient
                The gradient of the change in the objective or constraint with respect to each lambda.

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

        /// A vector of distances from the constraint manifold.
        std::vector<double> constraintDistances;

        /// A reference to a vector of optimum lambda values (to be found by the solver).
        std::vector<double>& lambdas;

        /// The effective time step.
        double& timeStep;

        /// The maximum displacement.
        double maxDisplacement;

        /// Whether to maximise the objective function.
        bool isMax;

        /// Whether the constraints are equality constraints.
        std::vector<bool> isEquality;

        /// The NLopt algorithm.
        nlopt::algorithm algorithm;

        /// A map between indices for active constraints.
        std::vector<unsigned int> indexMap;

        /// The boundary point displacement vector.
        std::vector<double> displacements;

        /// Negative lambda limits.
        std::vector<double> negativeLambdaLimits;

        /// Positive lambda limits.
        std::vector<double> positiveLambdaLimits;

        /// Scale factor for each function.
        std::vector<double> scaleFactors;

        /// Scaled constraint distances.
        std::vector<double> constraintDistancesScaled;

        /// Gradients of the objective and constraint functions.
        std::vector<std::vector<double> > gradients;

        /// Optimiser return code.
        nlopt::result returnCode;

        //! Compute scale factors.
        void computeScaleFactors();

        //! Compute constraint distances.
        /*! \param nCurrentConstraints
                The current number of active constraints.
         */
        void computeConstraintDistances(unsigned int);

        //! Compute lambda limits.
        void computeLambdaLimits();

        //! Compute the boundary point normal displacement.
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
        /*! \param gradient
                The gradient of the change in objective or constraint with respect to each lambda.

            \param index
                Function index, 0 = objective, 1, 2, 3, ... = constraints.
         */
        void computeGradients(std::vector<double>&, unsigned int index);

        //! Rescale displacements and lambda values if the CFL condition is violated.
        /*! \return
                The scale factor.
         */
        double rescaleDisplacements();
    };
}

#endif  /* _OPTIMISE_H */
