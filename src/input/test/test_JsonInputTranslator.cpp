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
#include <filesystem>

#include "JsonInputTranslator.h"

TEST(TestJsonInputTranslator, basic) {

    std::filesystem::path input_json = "../../../examples/unit_square/unit_square.json";

    JsonInputTranslator jit{input_json};

}

TEST(TestJsonInputTranslator, unit_square_multi_json) {


    std::filesystem::path input_json = "../../../examples/unit_square/unit_square.json";

    JsonInputTranslator jit{input_json};

    auto nu_map = jit.get_nu_map();

    for (auto [key, value] : nu_map){
        std::cout << key << "\t" << value->get_nu(1) << std::endl;
    }

}