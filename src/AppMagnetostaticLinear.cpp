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
#include "JsonInputTranslator.h"

using json = nlohmann::json;


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

    // Translate JSON
    JsonInputTranslator itrans{input, cli_source_map};

    auto nu_map = itrans.get_nu_map();
    auto f_map = itrans.get_f_map();
    auto dc_map = itrans.get_dc_map();
    auto per_map = itrans.get_per_map();
    auto postprocessors_cell = itrans.get_pp_cell();
    auto postprocessors_scalar = itrans.get_pp_scalar();
    auto mesh_path = itrans.get_mesh_filepath();
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

        // Create cell postprocessors
        std::map<std::string, ExpressionCellPostprocessor<2>*> user_expr_postprocessors;
        for (auto [key, val] : postprocessors_cell)
            user_expr_postprocessors[key] = new ExpressionCellPostprocessor<2>(val, nu_map, f_map);

        // Create scalar postprocessors
        std::map<std::string, ScalarPostprocessor<2>*> user_expr_sum_postprocessors;
        ScalarPostprocessorFactory<2> scalar_postprocessor_factory(nu_map, f_map);
        for (auto [key, val]: postprocessors_scalar)
            user_expr_sum_postprocessors[key] = scalar_postprocessor_factory.create(val);

        // Attach cell postprocessors to Exporters
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