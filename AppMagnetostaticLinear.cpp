//
// Created by gordan on 5/12/23.
//

#include <iostream>
#include <fstream>
#include <filesystem>
#include "LinearSolver.h"
#include "nlohmann/json.hpp"
#include <cxxopts.hpp>
#include <vector>

#include "ExportVtu.h"
#include "MagneticFluxPostprocessor.h"
#include "MatIDPostprocessor.h"
#include "MagneticEnergyDensityPostprocessor.h"
#include "MagneticEnergyPostprocessor.h"
#include "misc.h"
#include "exprtk.hpp"


using json = nlohmann::json;

int main(int argc, char* argv[]){

    cxxopts::Options options("aml", "Application for magnetostatic simulations.");

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
        std::cout << "No input file given.\nSee 'aml --help' for usage." << std::endl;
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
    json input_data = json::parse(ifs_input);
    std::string mesh_filepath{input_data.at("mesh_path")};

    std::cout << "Mesh file: " << mesh_filepath << std::endl;

    auto mesh_id_data = input_data.at("mesh_id");
    auto material_data = input_data.at("material");
    auto boundary_data = input_data.at("boundary");
    auto source_data = input_data.at("source");

    // modify source data if needed:
    // TODO: Rewrite possibly
    std::cout << "Sources: " << std::endl;
    for (auto& [key, val] : source_data.items()){
        if (cli_source_map.count(key))
            val = cli_source_map[key];
    }

    for (auto& [key, val] : source_data.items()){
        std::cout << key << ": " << val << std::endl;
    }

    std::unordered_map<int, std::variant<double, std::pair<double, double>>> f_map;
    std::unordered_map<int, double> nu_map;
    std::unordered_map<int, double> dc_map;

    static const double pi = 3.141592653589793238462643383279502;

    for (auto& mesh_el_data : mesh_id_data.items()){
        int mat_id{std::stoi(mesh_el_data.key())};
        if (mesh_el_data.value().contains("material")){
            nu_map.insert({mat_id, material_data.at(mesh_el_data.value().at("material")).at("nu")});
        }
        if (mesh_el_data.value().contains("source")){

            // if number take number; if string evaluate expression
            auto source_d = source_data.at(mesh_el_data.value().at("source"));
            double source_val = 0;
            if (source_d.is_number()){
                source_val = source_d;
            }
            else if (source_d.is_string()){
                std::string user_expr_string = source_d;
                std::cout << "user_expr_string: " << user_expr_string << std::endl;
                symbol_table_t symbol_table;
                expression_t expression;
                parser_t parser;

                symbol_table.add_constant("pi", pi);
                expression.register_symbol_table(symbol_table);
                parser.compile(user_expr_string, expression);
                source_val = expression.value();
                std::cout << std::setprecision(12);
                std::cout << "Evaluated " << mesh_el_data.value().at("source") << ": " << source_val << std::endl;
            }

            // If mesh data contains angle, the source is vector Hc. Magnitude in Material data and direction from mesh data.
            if (mesh_el_data.value().contains("angle")){
                double Hc = source_val;
                double angle = mesh_el_data.value().at("angle");
                double Hc_x = Hc*std::cos(angle);
                double Hc_y = Hc*std::sin(angle);

                std::pair<double, double> Hc_vec{Hc_x, Hc_y};
                f_map.insert({mat_id, Hc_vec});
            }

            // Otherwise, the source is J (current density)
            else {
                double J = source_val;
                f_map.insert({mat_id, J});
            }
        }

        // If there is no source, set f to 0.
        else{
            f_map.insert({mat_id, 0.0});
        }
        if (mesh_el_data.value().contains("boundary")){
            dc_map.insert({mat_id, boundary_data.at(mesh_el_data.value().at("boundary"))});
        }
    }

    // Initialize Solver and solver
    LinearSolver<2> solver;
    solver.read_mesh(mesh_filepath);
    solver.setup_system();
    solver.set_nu_map(nu_map);
    solver.set_f_map(f_map);
    solver.set_dc_map(dc_map);
    solver.assemble_system();
    solver.solve();

    // Visualization
    // Export
    MagneticFluxPostprocessor<2> bx_postprocessor(0, 0);
    MagneticFluxPostprocessor<2> by_postprocessor(0, 1);
    MagneticFluxPostprocessor<2> b_abs_postprocessor(0);
    MagneticEnergyDensityPostprocessor<2> e_density(nu_map);
    MagneticEnergyPostprocessor<2> e_cell(nu_map);
    MatIDPostprocessor<2> mat_id_postprocessor;

    std::vector<double> e_cells;
    e_cell.process(solver.get_triangulation(), solver.get_solution(), solver.get_fe(), e_cells);
    auto E_total = std::reduce(e_cells.begin(), e_cells.end());
    std::cout << "E total: " << E_total << std::endl;

    ExportVtu<2> export_vtu(solver.get_triangulation(), solver.get_rhs(), solver.get_solution(), solver.get_fe());
    export_vtu.attach_postprocessor(&mat_id_postprocessor, "MatID");
    export_vtu.attach_postprocessor(&b_abs_postprocessor, "|B| [T]");
    export_vtu.attach_postprocessor(&bx_postprocessor, "Bx [T]");
    export_vtu.attach_postprocessor(&by_postprocessor, "By [T]");
    export_vtu.attach_postprocessor(&e_density, "E [J/m3]");

    export_vtu.write(output);
    std::cout << "Output written to " << output.concat(".vtu") << std::endl;

    return 0;
}