//
// Created by gordan on 6/26/23.
//

#ifndef SOLVER_CONFIGPARSER_H
#define SOLVER_CONFIGPARSER_H

#include <iostream>
#include <fstream>
#include <filesystem>
#include "nlohmann/json.hpp"

using json = nlohmann::json;

class ConfigParser {

public:
    ConfigParser(std::filesystem::path input){
        root = input.parent_path();

        std::ifstream ifs_input;
        ifs_input.open(input);
        if(ifs_input.fail()){
            throw std::runtime_error("Failed to open input file.");
        }
        json input_data = json::parse(ifs_input);
        parse_json_data(input_data);
    }

    void parse_json_data(json data){
        for (auto [key, value] : data.items()){
            if (value.is_string()){
                std::string value_str = value;
                if (value_str.find(".json") != std::string::npos){
                    std::filesystem::path pathson{value_str};
                    if (pathson.is_relative()){
                        pathson = root / pathson;
                    }
                    std::ifstream ifs_input;
                    ifs_input.open(pathson);
                    if(ifs_input.fail()){
                        throw std::runtime_error("Failed to open " + key + " file: " + value_str);
                    }
                    parse_json_data(json::parse(ifs_input));
                }
                else
                    top_data[key] = value;
            }
            else
                top_data[key] = value;
        }
    }

    json get_top_data(){
        return top_data;
    }

    std::filesystem::path get_root(){
        return root;
    }


private:
    json top_data;
    std::filesystem::path root;
};


#endif //SOLVER_CONFIGPARSER_H
