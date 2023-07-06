//
// Created by gordan on 7/5/23.
//

#include "LinearBHCurve.h"

LinearBHCurve::LinearBHCurve(double nu) {
    this->nu = nu;
}

double LinearBHCurve::get_nu(double) {
    return this->nu;
}

double LinearBHCurve::get_nu_prime(double) {
    return 0;
}
