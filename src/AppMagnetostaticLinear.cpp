//
// Created by gordan on 5/12/23.
//

#include <iostream>
#include <fstream>
#include <filesystem>
#include <vector>

#include "nlohmann/json.hpp"
#include <cxxopts.hpp>
#include "exprtk.hpp"

#include "LinearSolver.h"
#include "export/ExportVtu.h"
#include "processors/ExpressionCellPostprocessor.h"
#include "processors/ScalarPostprocessorFactory.h"
#include "misc.h"
#include "ConfigParser.h"
#include "NuCurveFactory.h"
#include "FSourceFactory.h"

using json = nlohmann::json;

typedef std::unordered_map<int, NuCurve*> t_nu_map;
typedef std::unordered_map<int, std::variant<FSource*, std::pair<double, double>>> t_f_map;
typedef std::unordered_map<int, double> t_dc_map;
typedef std::unordered_map<std::string, std::vector<unsigned int>> t_per_map;
typedef std::unordered_map<std::string, double> t_cli_source_map;
typedef std::map<std::string, std::string> t_postprocessor_strings;

class InputTranslator{
    public:
        InputTranslator(json, t_cli_source_map);
        t_nu_map get_nu_map();
        t_f_map get_f_map();
        t_dc_map get_dc_map();
        t_per_map get_per_map();
        t_postprocessor_strings get_pp_cell();
        t_postprocessor_strings get_pp_scalar();

    private:
        json input_data;
        t_f_map f_map;
        t_nu_map nu_map;
        t_dc_map dc_map;
        t_per_map per_map;
        t_postprocessor_strings postprocessor_strings_cell;
        t_postprocessor_strings postprocessor_strings_scalar;

};

t_dc_map InputTranslator::get_dc_map() {
    return dc_map;
}

t_f_map InputTranslator::get_f_map() {
    return f_map;
}

t_nu_map InputTranslator::get_nu_map() {
    return nu_map;
}

t_per_map InputTranslator::get_per_map() {
    return per_map;
}

t_postprocessor_strings InputTranslator::get_pp_cell() {
    return postprocessor_strings_cell;
}

t_postprocessor_strings InputTranslator::get_pp_scalar() {
    return postprocessor_strings_scalar;
}

InputTranslator::InputTranslator(json input_data, t_cli_source_map cli_source_map) {

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

int main(int argc, char* argv[]){

    cxxopts::Options options("magas", "Application for magnetostatic simulations.");

    options.add_options()
            ("input", "Input config file (.json)", cxxopts::value<std::string>())
            ("o, output", "Visualization output", cxxopts::value<std::string>())
            ("m, mesh", "Mesh file (.gmsh)")
            ("s, sources", "Sources", cxxopts::value<std::vector<std::string>>())
            ("h, help", "Print usage")
    ;

    options.positional_help("file...");
    options.show_positional_help();
    options.parse_positional({"input"});
    options.allow_unrecognised_options();

    auto result = options.parse(argc, argv);

    if (result.count("help")){
        std::cout << options.help() << std::endl;
        exit(0);
    }

    if (!result.count("input")){
        std::cout << "No input file given.\nSee 'magas --help' for usage." << std::endl;
        exit(0);
    }

    std::filesystem::path input = result["input"].as<std::string>();
    std::filesystem::path output = input.filename().replace_extension();
    if (result.count("output")){
        output = result["output"].as<std::string>();
        output.replace_extension();
    }

    typedef exprtk::symbol_table<double> symbol_table_t;
    typedef exprtk::expression<double>   expression_t;
    typedef exprtk::parser<double>       parser_t;

    std::unordered_map<std::string, double> cli_source_map;
    if (result.count("sources")){
        for (const auto& src_param_str : result["sources"].as<std::vector<std::string>>()){
            std::vector<std::string> v = split(src_param_str, "=");
            std::string user_expr_string = v[1];
            symbol_table_t symbol_table;
            expression_t expression;
            parser_t parser;
            parser.compile(user_expr_string, expression);
            cli_source_map[v[0]] = expression.value();
        }
    }

    std::ifstream ifs_input;
    ifs_input.open(input);
        if(ifs_input.fail()){
            throw std::runtime_error("Failed to open input file.");
        }

    // Parse JSON
    ConfigParser cp{input};
    json input_data = cp.get_top_data();
    auto json_dir = cp.get_root();

    std::filesystem::path mesh_path{input_data.at("mesh_path")};
    if (mesh_path.is_relative()){
        mesh_path = json_dir / mesh_path;
    }

    InputTranslator itrans{input_data, cli_source_map};

    auto nu_map = itrans.get_nu_map();
    auto f_map = itrans.get_f_map();
    auto dc_map = itrans.get_dc_map();
    auto per_map = itrans.get_per_map();
    auto postprocessors_cell = itrans.get_pp_cell();
    auto postprocessors_scalar = itrans.get_pp_scalar();
    // Initialize Solver and solver
    try{
        LinearSolver<2> solver;
        std::cout << "Reading mesh..." << std::endl;
        solver.read_mesh(mesh_path);
        std::cout << "Setting up maps..." << std::endl;
        solver.set_nu_map(nu_map);
        solver.set_f_map(f_map);
        solver.set_dc_map(dc_map);
        solver.set_per_map(per_map);
        std::cout << "Setting up system..." << std::endl;
        solver.setup_system();
        std::cout << "Assembling system..." << std::endl;
        solver.assemble_system();
        std::cout << "Solving system..." << std::endl;
        solver.solve();

        // Create vector postprocessors
        std::map<std::string, ExpressionCellPostprocessor<2>*> user_expr_postprocessors;
        for (auto [key, val] : postprocessors_cell)
            user_expr_postprocessors[key] = new ExpressionCellPostprocessor<2>(val, nu_map, f_map);

        // Create scalar postprocessors
        std::map<std::string, ScalarPostprocessor<2>*> user_expr_sum_postprocessors;
        ScalarPostprocessorFactory<2> scalar_postprocessor_factory(nu_map, f_map);
        for (auto [key, val]: postprocessors_scalar)
            user_expr_sum_postprocessors[key] = scalar_postprocessor_factory.create(val);

        // Attach vector postprocessors to Exporters
        ExportVtu<2> export_vtu(solver.get_triangulation(), solver.get_rhs(), solver.get_solution(), solver.get_fe());
        for (auto [key, val] : user_expr_postprocessors)
            export_vtu.attach_postprocessor(val, key);

        // Perform postprocessing of scalar postprocessors
        std::unordered_map<std::string, double> results_map;
        double result_sum = 0;
        for (auto [key, val] : user_expr_sum_postprocessors){
            val->process(solver.get_triangulation(), solver.get_solution(), solver.get_fe(), result_sum);
            results_map[key] = result_sum;
            delete val;
        }

        // Perform vector postprocessing and export to vtu.
        // TODO: separate pefrom and write.
        export_vtu.write(output);
        std::cout << "Output written to " << output.concat(".vtu") << std::endl;
        for (auto [key, val] : results_map){
            std::cout << key << " = " << val << std::endl;
        }
        return 0;

    } catch (std::exception &exc){
        std::cerr << "Excpetion throw:" << std:: endl;
        std::cerr << exc.what() << std::endl;
        return 1;
    }
    catch (...){
        std::cerr << "Unknown exception! " << std::endl;
        return 1;
    }
}