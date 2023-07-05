//
// Created by gordan on 7/5/23.
//

#ifndef MAGAS_LINEARBHCURVE_H
#define MAGAS_LINEARBHCURVE_H

#include "../../solver/include/BHCurve.h"

class LinearBHCurve : public BHCurve {

public:
    LinearBHCurve(double nu);
    double get_nu(double) override;
    double get_nu_prime(double) override;

private:
    double nu;
};


#endif //MAGAS_LINEARBHCURVE_H
