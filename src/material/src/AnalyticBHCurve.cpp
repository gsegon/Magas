//
// Created by gordan on 7/5/23.
//

#include "../include/AnalyticNuCurve.h"
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
