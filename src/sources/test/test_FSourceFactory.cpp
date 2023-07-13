//
// Created by gordan on 3/12/23.
//
#include <gtest/gtest.h>

#include "FSourceFactory.h"

TEST(TestFSourceFactory, fconst) {

    FSourceFactory ff;

    auto f_const = ff.create(5);
    ASSERT_NEAR(f_const->get_value(0, 0), 5, 1e-8);

}

TEST(TestFSourceFactory, fexpr) {

    FSourceFactory ff;

    auto fx = ff.create("x");
    ASSERT_NEAR(fx->get_value(5, 4), 5, 1e-8);

    auto fx2 = ff.create("x^2");
    ASSERT_NEAR(fx2->get_value(5, 4), 25, 1e-8);

    auto fy = ff.create("y");
    ASSERT_NEAR(fy->get_value(5, 4), 4, 1e-8);

    auto fy2 = ff.create("y^2");
    ASSERT_NEAR(fy2->get_value(5, 4), 16, 1e-8);

}
