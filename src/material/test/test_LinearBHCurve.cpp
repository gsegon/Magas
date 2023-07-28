//
// Created by gordan on 3/12/23.
//
#include <gtest/gtest.h>

#include "LinearNuCurve.h"

TEST(LinearBHCurve, basic) {

    LinearNuCurve bh1{795774.715025};

    ASSERT_EQ(bh1.get_nu(0.5), 795774.715025);
    ASSERT_EQ(bh1.get_nu_prime(0.5), 0.0);
}

TEST(LinearBHCurve, basic2) {

    LinearNuCurve bh1{100};

    ASSERT_EQ(bh1.get_nu(0.5), 100);
    ASSERT_EQ(bh1.get_nu_prime(0.5), 0.0);
}

