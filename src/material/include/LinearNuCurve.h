//
// Created by gordan on 7/5/23.
//

#ifndef MAGAS_LINEARNUCURVE_H
#define MAGAS_LINEARNUCURVE_H

#include "NuCurve.h"

class LinearNuCurve : public NuCurve {

public:
    LinearNuCurve(double nu);
    double get_nu(double) override;
    double get_nu_prime(double) override;
    double get_energy(double) override;
    double get_coenergy(double) override;

private:
    double nu;
};


#endif //MAGAS_LINEARNUCURVE_H
