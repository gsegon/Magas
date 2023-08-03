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

#include <gtest/gtest.h>

#include "AnalyticNuCurve.h"

TEST(AnalyticBHCurve, basic) {

    AnalyticNuCurve bhanalytic1;

    std::cout << bhanalytic1.get_nu(2.0) << std::endl;
    std::cout << bhanalytic1.get_nu_prime(2.0) << std::endl;

    ASSERT_NEAR(bhanalytic1.get_nu(2.0), 21744.9, 1e-1);
    ASSERT_NEAR(bhanalytic1.get_nu_prime(2.0), 94210, 1e-1);
}
