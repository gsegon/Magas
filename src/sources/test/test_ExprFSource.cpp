//
// Created by gordan on 3/12/23.
//
#include <gtest/gtest.h>

#include "ExprFSource.h"

TEST(ExprFSource, basic) {

    ExprFSource f{"x"};
    ASSERT_NEAR(f.get_value(5, 4), 5, 1e-8);

    ExprFSource f2{"x^2"};
    ASSERT_NEAR(f2.get_value(5, 4), 25, 1e-8);

    ExprFSource f3{"y"};
    ASSERT_NEAR(f3.get_value(5, 4), 4, 1e-8);

    ExprFSource f4{"y^2"};
    ASSERT_NEAR(f4.get_value(5, 4), 16, 1e-8);

}
