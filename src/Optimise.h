/*
 * File:    Optimise.h
 * Author:  lester
 */

#ifndef _OPTIMISE_H
#define _OPTIMISE_H

#include <nlopt.hpp>

#include "Boundary.h"

// ASSOCIATED DATA TYPES

struct NLOptWrapper
{
    unsigned int index;     //!< Function index (0 = objective, 1, 2, ... = constraints).
    void* callback;         //!< Pointer to callback function wrapper.
};

//! NLOpt callback function wrapper.
/*! \param lambda
        The current lambda value for each function.

    \param gradient
        The gradient of the change in objective or constraint with respect to each lambda.

    \param data
        Void pointer to NLOptWrapper data.
 */
double callbackWrapper(const std::vector<double>& lambda, std::vector<double>& gradient, void* data);

/*!\brief A class to solve for the optimum boundary movement vector (velocity).
 */
class Optimise
{
public:
    //! Constructor.
    /*! \param boundaryPoints_
            A reference to a vector of boundary points.

        \param constraintDistances_
            Distance from violation for each constraint.

        \param lambdas_
            The optimum lambda values. This array is modified.

        \param velocities_
            The optimum boundary movement vector (velocity). This array is modified.
     */
    Optimise(const std::vector<BoundaryPoint>&, const std::vector<double>&,
        std::vector<double>&, std::vector<double>&);

    //! NLOpt callback function.
    /*! \param lambda
            The current lambda value for each function.

        \param gradient
            The gradient of the change in objective or constraint with respect to each lambda.

        \param index
            Function index, 0 = objective, 1, 2, 3, ... = constraints.
     */
    double callback(const std::vector<double>&, std::vector<double>&, unsigned int);

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
};

#endif  /* _OPTIMISE_H */
