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

    ExprFSource f5{"5.35"};
    ASSERT_NEAR(f5.get_value(5, 4), 5.35, 1e-8);

    ExprFSource f6{"0*5.35"};
    ASSERT_NEAR(f6.get_value(5, 4), 0, 1e-8);

}
