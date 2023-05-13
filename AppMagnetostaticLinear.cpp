//
// Created by gordan on 5/12/23.
//

#include <iostream>
#include <fstream>
#include "LinearSolver.h"
#include "nlohmann/json.hpp"

using json = nlohmann::json;

int main(int argc, char* argv[]){

    std::ifstream input;

    if (argc <= 1){
        std::cout << "Error: No input file given." << std::endl;
        return -1;
    }
    if (argc > 1){
        input.open(argv[1]);
        if(input.fail()){
            throw std::runtime_error("Failed to open input file");
        }
        std::cout << "Input file " << argv[1] << " loaded." << std::endl;
    }

    // Parse JSON
    json input_data = json::parse(input);

    std::string mesh_filepath{input_data.at("mesh_path")};
    auto mesh_id_data = input_data.at("mesh_id");
    auto material_data = input_data.at("material");
    auto boundary_data = input_data.at("boundary");
    auto source_data = input_data.at("source");

    std::unordered_map<int, double> f_map;
    std::unordered_map<int, double> nu_map;
    std::unordered_map<int, double> dc_map;

    for (auto& mesh_data : mesh_id_data.items()){
        int mat_id{std::stoi(mesh_data.key())};
        if (mesh_data.value().contains("material")){
            nu_map.insert({mat_id, material_data.at(mesh_data.value().at("material")).at("nu")});
        }
        if (mesh_data.value().contains("source")){
            f_map.insert({mat_id, source_data.at(mesh_data.value().at("source"))});
        }
        if (mesh_data.value().contains("boundary")){
            dc_map.insert({mat_id, boundary_data.at(mesh_data.value().at("boundary"))});
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

    // Output results
    solver.output_results("alm_output");

    return 0;
}