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


#include "JsonInputTranslator.h"
#include "ConfigParser.h"
#include "InterpolatedNuCurve.h"


t_dc_map JsonInputTranslator::get_dc_map() {
    return dc_map;
}

t_f_map JsonInputTranslator::get_f_map() {
    return f_map;
}

t_nu_map JsonInputTranslator::get_nu_map() {
    return nu_map;
}

t_per_map JsonInputTranslator::get_per_map() {
    return per_map;
}

t_rot_map JsonInputTranslator::get_rot_map() {
    return rot_map;
}

t_postprocessor_strings JsonInputTranslator::get_pp_cell() {
    return postprocessor_strings_cell;
}

t_postprocessor_strings JsonInputTranslator::get_pp_scalar() {
    return postprocessor_strings_scalar;
}

std::filesystem::path JsonInputTranslator::get_mesh_filepath() {
    return mesh_filepath;
}

JsonInputTranslator::JsonInputTranslator(std::filesystem::path input) {

    ConfigParser cp{input};
    json input_data = cp.get_top_data();
    auto json_dir = cp.get_root();

    mesh_filepath = (std::filesystem::path)input_data.at("mesh_path");
    if (mesh_filepath.is_relative()){
        mesh_filepath = json_dir / mesh_filepath;
    }

    auto material_data = input_data.at("material");
    auto boundary_data = input_data.at("boundary");
    auto source_data = input_data.at("source");
    json rotation_data;
    if (input_data.contains("rotation")){
        rotation_data = input_data.at("rotation");
    }

    auto boundary_id_data = input_data.at("boundary_id");
    auto mesh_id_data = input_data.at("mesh_id");
    auto postprocess_data = input_data.at("postprocess");
    auto postprocess_sum_data = input_data.at("postprocess_sum");

    for (auto [key, val] : postprocess_data.items())
        postprocessor_strings_cell[key] = val;

    for (auto [key, val] : postprocess_sum_data.items())
        postprocessor_strings_scalar[key] = val;

    // Add boundary values to dc_map
    for (auto& boundary_el_data : boundary_id_data.items()) {
        int boundary_id{std::stoi(boundary_el_data.key())};
        auto boundary_value = boundary_data.at(boundary_el_data.value().at("boundary"));
        if (boundary_value.is_number())
            dc_map.insert({boundary_id, boundary_value});
        if (boundary_value.is_string())
            per_map[boundary_value].push_back(boundary_id);
    }

    // Add material coefficients to 'nu_map' and sources to 'f_map'
    NuCurveFactory bh_factory;
    FSourceFactory fsf;
    std::map<std::string, std::vector<unsigned int>> rot_to_mat_ids;
    for (auto& mesh_el_data : mesh_id_data.items()){
        int mat_id{std::stoi(mesh_el_data.key())};
        if (mesh_el_data.value().contains("material")){
            auto value1 = material_data.at(mesh_el_data.value().at("material")).at("nu");
            if (value1.is_number()) nu_map.insert({mat_id, bh_factory.create((double)value1)});
            if (value1.is_string()){
                if (((std::string)value1).find(".csv") != std::string::npos){
                    filesystem::path bhpath{value1};
                    if (bhpath.is_relative()){
                        bhpath = json_dir / bhpath;
                    }
                    nu_map.insert({mat_id, bh_factory.create(bhpath)});
                }
                else{
                    nu_map.insert({mat_id, bh_factory.create((string)value1)});
                }
            }
        }

        if (mesh_el_data.value().contains("rotation")) {
            auto value1 =mesh_el_data.value().at("rotation");
            rot_to_mat_ids[value1].push_back(mat_id);
//            auto offset = rotation_data.at(value1);
        }

        if (mesh_el_data.value().contains("source")){

            // if number take number; if string evaluate expression
            auto source_d = source_data.at(mesh_el_data.value().at("source"));
            double source_val = 0;
            if (!mesh_el_data.value().contains("angle")){

                if (source_d.is_number()){
                    source_val = source_d;
                    f_map.insert({mat_id, fsf.create(source_val)});
                }
                else if (source_d.is_string()){
                    std::string user_expr_string = source_d;
                    f_map.insert({mat_id, fsf.create(user_expr_string)});
                }
            }

            // If mesh data contains angle, the source is vector Hc. Magnitude in Material data and direction from mesh data.
            if (mesh_el_data.value().contains("angle")){
                if (source_d.is_number()){
                    source_val = source_d;
                }
                double Hc = source_val;
                double angle = mesh_el_data.value().at("angle");
                double Hc_x = Hc*std::cos(angle);
                double Hc_y = Hc*std::sin(angle);

                std::pair<double, double> Hc_vec{Hc_x, Hc_y};
                f_map.insert({mat_id, Hc_vec});
            }
        }
            // If there is no source, set f to 0.
        else{
            f_map.insert({mat_id, fsf.create(0)});
        }
    }

    if (input_data.contains("rotation")) {
        for (auto &rot_data: rotation_data.items()) {
            std::pair<unsigned int, unsigned int> rot_pair{rot_to_mat_ids[rot_data.key()][0],
                                                           rot_to_mat_ids[rot_data.key()][1]};
            rot_map.insert({rot_pair, rot_data.value()});
        }
    }
}

bool JsonInputTranslator::is_nonlinear() {
    for (auto [key, val] : nu_map){
        if (typeid(*val) == typeid(InterpolatedNuCurve)){
            nonlinear = true;
        }
    }
    return nonlinear;
}