//
// Created by gordan on 3/12/23.
//
#include <gtest/gtest.h>

#include "../include/BHCurveFactory.h"

TEST(TestBHCurveFactory, basic) {


    BHCurveFactory bhcf;
    auto a = bhcf.create(10);
    a->get_nu(0.5);

    ASSERT_EQ( a->get_nu(0.5), 10);

}
