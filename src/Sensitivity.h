/*
 * File:	Sensitivity.h
 * Author:	lester
 */

#ifndef _SENSITIVITY_H
#define _SENSITIVITY_H

#include <functional>

#include "Boundary.h"

/*! \file Sensitivity.h
    \brief A class for calculating finite-difference boundary point sensitivities.
 */

namespace lsm
{
    //! Calculate the value of a function for a small displacement of a boundary point.
    /*! \param point
            A reference to the boundary point.

        \return
            The value of the function.
     */
    typedef std::function<double (const BoundaryPoint&)> SensitivityCallback;

    /*! \brief A class for calculating finite-difference boundary point sensitivities.
    */
    class Sensitivity
    {
    public:
        //! Constructor.
        /*! \param delta_
                The finite difference derivative length.
        */
        Sensitivity(double delta_ = 0.1);

        //! Compute finite difference sensitivity for an arbitrary function.
        /*! \param point
                A reference to the boundary point.

            \param callback
                A reference to the sensitivity callback function.

            \return
                The finite-difference sensitivity.
         */
        double computeSensitivity(BoundaryPoint&, SensitivityCallback&);

    private:
        /// The finite-difference derivative length (in units of the grid spacing).
        double delta;
    };
}

#endif	/* _SENSITIVITY_H */
