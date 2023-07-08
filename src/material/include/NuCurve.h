//
// Created by gordan on 7/5/23.
//

#ifndef MAGAS_NUCURVE_H
#define MAGAS_NUCURVE_H

class NuCurve{

public:
    virtual double get_nu(double) = 0;
    virtual double get_nu_prime(double) = 0;
    virtual double get_energy(double) = 0;
    virtual double get_coenergy(double) = 0;
};

#endif //MAGAS_NUCURVE_H
