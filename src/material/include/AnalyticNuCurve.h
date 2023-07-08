//
// Created by gordan on 7/5/23.
//

#ifndef MAGAS_ANALYTICNUCURVE_H
#define MAGAS_ANALYTICNUCURVE_H

#include "NuCurve.h"

class AnalyticNuCurve : public NuCurve {

public:
    AnalyticNuCurve();
    double get_nu(double) override;
    double get_nu_prime(double) override;
    double get_energy(double) override;
    double get_coenergy(double) override;

private:
    double nu_0;
    double alpha;
    double beta;
    double theta;
    double b_sat;
};


#endif //MAGAS_ANALYTICNUCURVE_H
