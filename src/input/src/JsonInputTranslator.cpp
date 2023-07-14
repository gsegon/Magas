//
// Created by gordan on 7/14/23.
//

#include "JsonInputTranslator.h"
#include "ConfigParser.h"


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

t_postprocessor_strings JsonInputTranslator::get_pp_cell() {
    return postprocessor_strings_cell;
}

t_postprocessor_strings JsonInputTranslator::get_pp_scalar() {
    return postprocessor_strings_scalar;
}

std::filesystem::path JsonInputTranslator::get_mesh_filepath() {
    return mesh_filepath;
}

JsonInputTranslator::JsonInputTranslator(std::filesystem::path input, t_cli_source_map cli_source_map) {

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
    auto boundary_id_data = input_data.at("boundary_id");
    auto mesh_id_data = input_data.at("mesh_id");

    auto postprocess_data = input_data.at("postprocess");
    auto postprocess_sum_data = input_data.at("postprocess_sum");

    for (auto [key, val] : postprocess_data.items())
        postprocessor_strings_cell[key] = val;

    for (auto [key, val] : postprocess_sum_data.items())
        postprocessor_strings_scalar[key] = val;

    // modify source data if needed:
    // TODO: Rewrite possibly
    for (auto& [key, val] : source_data.items()){
        if (cli_source_map.count(key))
            val = cli_source_map[key];
    }

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
    NuCurveFactory bhcf;
    FSourceFactory fsf;
    for (auto& mesh_el_data : mesh_id_data.items()){
        int mat_id{std::stoi(mesh_el_data.key())};
        if (mesh_el_data.value().contains("material")){
            auto value1 = material_data.at(mesh_el_data.value().at("material")).at("nu");
            if (value1.is_number()) nu_map.insert({mat_id, bhcf.create((double)value1)});
            if (value1.is_string()) nu_map.insert({mat_id, bhcf.create((string)value1)});
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
}