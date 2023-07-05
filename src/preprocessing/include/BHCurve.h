//
// Created by gordan on 7/5/23.
//

#ifndef MAGAS_BHCURVE_H
#define MAGAS_BHCURVE_H

class BHCurve{

public:
    virtual double get_nu(double) = 0;
    virtual double get_nu_prime(double) = 0;
};

#endif //MAGAS_BHCURVE_H
