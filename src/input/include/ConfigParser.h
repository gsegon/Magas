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
