//
// Created by gordan on 3/12/23.
//
#include <gtest/gtest.h>

#include "ConstFSource.h"

TEST(TestConstFSource, basic) {

    ConstFSource f{1};

    ASSERT_NEAR(f.get_value(0, 0), 1, 1e-8);

}
