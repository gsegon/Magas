//
// Created by gordan on 7/5/23.
//

#ifndef MAGAS_INTERPOLATEDNUCURVE_H
#define MAGAS_INTERPOLATEDNUCURVE_H

#include <vector>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "NuCurve.h"

class InterpolatedNuCurve : public NuCurve {

public:
    InterpolatedNuCurve(std::vector<double>, std::vector<double>);
    ~InterpolatedNuCurve();
    double get_nu(double) override;
    double get_nu_prime(double) override;

private:
    std::vector<double> b;
    std::vector<double> h;

    gsl_interp_accel *acc = nullptr;
    gsl_spline *spline = nullptr;

};


#endif //MAGAS_INTERPOLATEDNUCURVE_H
