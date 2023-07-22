//
// Created by gordan on 7/22/23.
//
//
// Created by gordan on 3/12/23.
//
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