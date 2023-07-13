//
// Created by gordan on 3/12/23.
//
#include <gtest/gtest.h>
#include <filesystem>

#include "ConfigParser.h"

TEST(TestConfigParser, basic) {

    std::filesystem::path input_json = "../../../examples/unit_square/unit_square.json";

    ConfigParser cp{input_json};

}

TEST(TestConfigParser, unit_square_multi_json) {


    std::filesystem::path input_json = "../../../examples/unit_square/unit_square.json";

    ConfigParser cp{input_json};

    auto data = cp.get_top_data();
    std::cout << "data: " << data << std::endl;


}

