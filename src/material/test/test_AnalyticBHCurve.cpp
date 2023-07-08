//
// Created by gordan on 3/12/23.
//
#include <gtest/gtest.h>

#include "AnalyticNuCurve.h"

TEST(AnalyticBHCurve, basic) {

    AnalyticNuCurve bhanalytic1;

    std::cout << bhanalytic1.get_nu(2.0) << std::endl;
    std::cout << bhanalytic1.get_nu_prime(2.0) << std::endl;

    ASSERT_NEAR(bhanalytic1.get_nu(2.0), 21744.9, 1e-1);
    ASSERT_NEAR(bhanalytic1.get_nu_prime(2.0), 94210, 1e-1);
}
