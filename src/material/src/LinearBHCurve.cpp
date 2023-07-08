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

double LinearNuCurve::get_coenergy(double b) {
    return get_energy(b);
}

double LinearNuCurve::get_energy(double b) {
    return nu*b*b/2;
}
