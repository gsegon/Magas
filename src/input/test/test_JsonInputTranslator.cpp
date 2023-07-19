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

TEST(TestJsonInputTranslator, rotation_benchmark_kelvin_1_json) {

    std::filesystem::path input_json = "/home/gordan/Programs/Magas/examples/rotation_benchmark_kelvin_1/rotation_benchmark_kelvin_1.json";

    JsonInputTranslator jit{input_json};

    auto nu_map = jit.get_nu_map();

    for (auto [key, value] : nu_map){
        std::cout << key << "\t" << value->get_nu(1) << std::endl;
    }

    for (auto [key, value] : jit.get_rot_map()){
        std::cout << key.first << "->" << key.second << "\toffset: " << value << std::endl;
    }

}

