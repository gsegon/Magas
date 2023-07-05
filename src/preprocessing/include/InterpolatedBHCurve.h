//
// Created by gordan on 7/5/23.
//

#ifndef MAGAS_INTERPOLATEDBHCURVE_H
#define MAGAS_INTERPOLATEDBHCURVE_H

#include <vector>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "../../solver/include/BHCurve.h"

class InterpolatedBHCurve : public BHCurve {

public:
    InterpolatedBHCurve(std::vector<double>, std::vector<double>);
    ~InterpolatedBHCurve();
    double get_nu(double) override;
    double get_nu_prime(double) override;

private:
    std::vector<double> b;
    std::vector<double> h;

    gsl_interp_accel *acc = nullptr;
    gsl_spline *spline = nullptr;

};


#endif //MAGAS_INTERPOLATEDBHCURVE_H
