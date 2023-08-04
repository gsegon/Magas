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
typedef std::map<std::pair<unsigned int, unsigned int>, int> t_rot_map;

class JsonInputTranslator{
public:
    JsonInputTranslator(std::filesystem::path);
    t_nu_map get_nu_map();
    t_f_map get_f_map();
    t_dc_map get_dc_map();
    t_per_map get_per_map();
    t_rot_map get_rot_map();
    t_postprocessor_strings get_pp_cell();
    t_postprocessor_strings get_pp_scalar();
    std::filesystem::path get_mesh_filepath();
    bool is_nonlinear();

private:
    json input_data;
    t_f_map f_map;
    t_nu_map nu_map;
    t_dc_map dc_map;
    t_per_map per_map;
    t_rot_map rot_map;
    t_postprocessor_strings postprocessor_strings_cell;
    t_postprocessor_strings postprocessor_strings_scalar;
    std::filesystem::path mesh_filepath;
    bool nonlinear=false;

};



#endif //MAGAS_JSONINPUTTRANSLATOR_H