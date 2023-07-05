//
// Created by gordan on 3/12/23.
//
#include <gtest/gtest.h>

#include "InterpolatedBHCurve.h"

TEST(InterpolatedBHCurve, basic) {


    std::vector<double> b{0, 0.1, 0.2, 0.3};
    std::vector<double> h{0, 1, 4, 16};

    InterpolatedBHCurve ibh(b, h);

    std::cout << "nu(0.15): " << ibh.get_nu(0.15) << std::endl;
    std::cout << "nu_prime(0.15): " << ibh.get_nu_prime(0.15) << std::endl;

}
