//
// Created by gordan on 7/5/23.
//

#ifndef MAGAS_ANALYTICBHCURVE_H
#define MAGAS_ANALYTICBHCURVE_H

#include "BHCurve.h"

class AnalyticBHCurve : public BHCurve {

public:
    AnalyticBHCurve();
    double get_nu(double) override;
    double get_nu_prime(double) override;

private:
    double nu_0;
    double alpha;
    double beta;
    double theta;
    double b_sat;
};


#endif //MAGAS_ANALYTICBHCURVE_H