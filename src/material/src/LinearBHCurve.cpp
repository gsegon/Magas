//
// Created by gordan on 7/5/23.
//

#include "LinearNuCurve.h"

LinearNuCurve::LinearNuCurve(double nu) {
    this->nu = nu;
}

double LinearNuCurve::get_nu(double) {
    return this->nu;
}

double LinearNuCurve::get_nu_prime(double) {
    return 0;
}
