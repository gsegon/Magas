//
// Created by gordan on 3/12/23.
//
#include <gtest/gtest.h>

#include "ConfigParser.h"

TEST(TestConfigParser, basic) {

    std::filesystem::path home = std::getenv("HOME");
    std::filesystem::path json = "Programs/solver/examples/unit_square/unit_square.json";
    std::string input_json{home/json};

    ConfigParser cp{input_json};

}

TEST(TestConfigParser, unit_square_multi_json) {

    std::filesystem::path home = std::getenv("HOME");
    std::filesystem::path json = "Programs/solver/examples/unit_square/unit_square.json";
    std::string input_json{home/json};

    ConfigParser cp{input_json};

    auto data = cp.get_top_data();
    std::cout << "data: " << data << std::endl;



}

