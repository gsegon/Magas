// Magas - Magnetostatic Analysis Suite
// Copyright (C) 2023  Gordan Segon
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
// USA

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
