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

#include "Sensitivity.h"

/*! \file Sensitivity.cpp
    \brief A class for calculating finite-difference boundary point sensitivities.
 */

namespace slsm
{
    Sensitivity::Sensitivity(double delta_) : delta(delta_)
    {
    }

    double Sensitivity::computeSensitivity(BoundaryPoint& point, SensitivityCallback& callback)
    {
        // Store the initial boundary point coordinates.
        Coord coord = point.coord;

        // Displace the boundary point in the positive normal direction.
        point.coord.x = coord.x + delta*point.normal.x;
        point.coord.y = coord.y + delta*point.normal.y;

        // Compute the new value of the function.
        double f1 = callback(point);

        // Displace the boundary point in the negative normal direction.
        point.coord.x = coord.x - delta*point.normal.x;
        point.coord.y = coord.y - delta*point.normal.y;

        // Compute the new value of the function.
        double f2 = callback(point);

        // Compute the finite-difference derivative.
        double sens = (f1 - f2) / (2.0 * delta);

        // Divide by boundary point length (sensitivity per unit length).
        sens /= point.length;

        // Reset boundary point coordinates.
        point.coord = coord;

        // Return sensitivity.
        return sens;
    }

    void Sensitivity::itoCorrection(Boundary& boundary, double temperature)
    {
        if (temperature == 0) return;

        // Compute the boundary normal vector.
        boundary.computeNormalVectors();

        using namespace std::placeholders;
        SensitivityCallback callback = std::bind(&Boundary::computePerimeter, boundary, _1);

        // Apply the deterministic Ito correction.
        for (unsigned int i=0;i<boundary.points.size();i++)
        {
            // Compute the local boundary point curvature.
            double curvature = computeSensitivity(boundary.points[i], callback);

            // Correct the objective sensitivity.
            boundary.points[i].sensitivities[0] -= (temperature * curvature) / (2.0 * boundary.points[i].length);
        }
    }
}
