// Magas - Magnetostatic Analysis Suite
// Copyright (C) 2023  Gordan Segon
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
// USA

#include "AnalyticNuCurve.h"
#include <cmath>

AnalyticNuCurve::AnalyticNuCurve() {
    nu_0 = 795774.715025;
    alpha = 2077.205761389225;
    beta = 5.289952851132246;
    b_sat = 2.2530727020352703;
}

namespace Analytic{
    double h_fun_1(double b){
        double alpha = 2077.205761389225;
        double beta = 5.289952851132246;

        return alpha*b+std::exp(beta*b)-1;
    }

    double h_fun(double b){
        double nu_0 = 795774.715025;
        double b_sat = 2.2530727020352703;
        double theta = -1638220.518181392;

        if (b < b_sat)
            return h_fun_1(b);
        else
            return nu_0*b + theta;
    }
}

double AnalyticNuCurve::get_nu(double b) {
    double alpha = 2077.205761389225;
    if (b == 0)
        return alpha;
    else
        return Analytic::h_fun(b)/b;
}

double AnalyticNuCurve::get_nu_prime(double b) {
    double nu_0 = 795774.715025;
    double b_sat = 2.2530727020352703;
    double theta = -1638220.518181392;
    double beta = 5.289952851132246;

    if (b == 0){
        return 0;
    }
    else if (b < b_sat)
        return beta*std::exp(beta*b)/b +(1-std::exp(beta*b))/(b*b);
    else
        return -theta/(b*b);
}

double AnalyticNuCurve::get_coenergy(double b) {
    return 0;
}

double AnalyticNuCurve::get_energy(double b) {
    return 0;
}

