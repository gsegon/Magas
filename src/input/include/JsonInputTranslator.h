//
// Created by gordan on 7/14/23.
//

#ifndef MAGAS_JSONINPUTTRANSLATOR_H
#define MAGAS_JSONINPUTTRANSLATOR_H

#include "nlohmann/json.hpp"
#include "NuCurveFactory.h"
#include "FSourceFactory.h"

using json = nlohmann::json;

typedef std::unordered_map<int, NuCurve*> t_nu_map;
typedef std::unordered_map<int, std::variant<FSource*, std::pair<double, double>>> t_f_map;
typedef std::unordered_map<int, double> t_dc_map;
typedef std::unordered_map<std::string, std::vector<unsigned int>> t_per_map;
typedef std::map<std::string, std::string> t_postprocessor_strings;

class JsonInputTranslator{
public:
    JsonInputTranslator(std::filesystem::path);
    t_nu_map get_nu_map();
    t_f_map get_f_map();
    t_dc_map get_dc_map();
    t_per_map get_per_map();
    t_postprocessor_strings get_pp_cell();
    t_postprocessor_strings get_pp_scalar();
    std::filesystem::path get_mesh_filepath();

private:
    json input_data;
    t_f_map f_map;
    t_nu_map nu_map;
    t_dc_map dc_map;
    t_per_map per_map;
    t_postprocessor_strings postprocessor_strings_cell;
    t_postprocessor_strings postprocessor_strings_scalar;
    std::filesystem::path mesh_filepath;

};



#endif //MAGAS_JSONINPUTTRANSLATOR_H
