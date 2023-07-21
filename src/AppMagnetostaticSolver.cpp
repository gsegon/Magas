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

#include "JsonInputTranslator.h"
#include "LinearSolver.h"
#include "export/ExportVtu.h"
#include "processors/ExpressionCellPostprocessor.h"
#include "processors/ScalarPostprocessorFactory.h"
#include "Solver.h"
#include "SolverFactory.h"

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

    // Translate JSON
    JsonInputTranslator itrans{input};

    // Translate maps
    auto nu_map = itrans.get_nu_map();
    auto f_map = itrans.get_f_map();
    auto dc_map = itrans.get_dc_map();
    auto per_map = itrans.get_per_map();
    auto rot_map = itrans.get_rot_map();
    auto postprocessors_cell = itrans.get_pp_cell();
    auto postprocessors_scalar = itrans.get_pp_scalar();
    auto mesh_path = itrans.get_mesh_filepath();

    // Initialize Solver and solver
    SolverFactory<2> sf;
    try{
        Solver<2>* solver = sf.create(itrans.is_nonlinear());
        std::cout << "::Initializing solver::" << std::endl;
        std::cout << "\tReading mesh...";
        solver->read_mesh(mesh_path);
        std::cout << "Done!" << std::endl;
        std::cout << "\tSetting up maps...";
        solver->set_nu_map(nu_map);
        solver->set_f_map(f_map);
        solver->set_dc_map(dc_map);
        solver->set_per_map(per_map);
        solver->set_rot_map(rot_map);
        std::cout << "Done!" << std::endl;
        solver->run();

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
        ExportVtu<2> export_vtu(solver->get_triangulation(), solver->get_rhs(), solver->get_solution(), solver->get_fe());
        for (auto [key, val] : user_expr_postprocessors)
            export_vtu.attach_postprocessor(val, key);

        // Perform postprocessing of scalar postprocessors
        std::unordered_map<std::string, double> results_map;
        double result_sum = 0;
        for (auto [key, val] : user_expr_sum_postprocessors){
            val->process(solver->get_triangulation(), solver->get_solution(), solver->get_fe(), result_sum);
            results_map[key] = result_sum;
        }

        json results;
        std::cout << "::Postprocessing results::" << std::endl;
        for (auto [key, val] : results_map){
            std::cout << "\t" << key << " = " << val << std::endl;
            results[key] = val;
        }

        std::ofstream o("results- " + (std::string)output + ".json");
        o << std::setw(4) << results << std::endl;

        // Perform vector postprocessing and export to vtu.
        // TODO: separate perform and write.
        export_vtu.write(output);
        std::cout << "::Export results::" << std::endl;
        std::cout << "\tOutput written to " << output.concat(".vtu").string();
        std::cout << " (" << (std::filesystem::current_path() / output)<< ")" << std::endl;

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