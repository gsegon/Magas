//
// Created by gordan on 5/12/23.
//

#include <iostream>
#include <fstream>
#include "PicardSolver.h"
#include "nlohmann/json.hpp"
#include "ExportVtu.h"
#include "MagneticFluxPostprocessor.h"
#include "MatIDPostprocessor.h"

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

    std::cout << "Mesh file: " << mesh_filepath << std::endl;

    auto mesh_id_data = input_data.at("mesh_id");
    auto material_data = input_data.at("material");
    auto boundary_data = input_data.at("boundary");
    auto source_data = input_data.at("source");

    std::unordered_map<int, std::any> nu_map;
    std::unordered_map<int, double> f_map;
    std::unordered_map<int, double> dc_map;

    for (auto& mesh_data : mesh_id_data.items()){
        int mat_id{std::stoi(mesh_data.key())};
        if (mesh_data.value().contains("material")){
            auto value = material_data.at(mesh_data.value().at("material")).at("nu");

            if (value.type() == json::value_t::string)
                nu_map.insert({mat_id, value});
            else
                nu_map.insert({mat_id, (double) value});
        }
        if (mesh_data.value().contains("source")){
            f_map.insert({mat_id, source_data.at(mesh_data.value().at("source"))});
        }
        else{
            f_map.insert({mat_id, 0.0});
        }
        if (mesh_data.value().contains("boundary")){
            dc_map.insert({mat_id, boundary_data.at(mesh_data.value().at("boundary"))});
        }
    }

    // Initialize Solver and solver
    constexpr double mu_0 = 1.2566370621219e-6;
    constexpr double nu_core = 1/(2500*mu_0);

    PicardSolver<2> solver;
    solver.read_mesh(mesh_filepath);
    solver.setup_cell_nu_history();
    solver.setup_system();
    solver.set_nu_map(nu_map);
    solver.set_f_map(f_map);
    solver.set_dc_map(dc_map);
    solver.initialize_cell_nu_history(nu_core);
    solver.solve_nonlinear(10);

    // Visualization
    // Export
    MagneticFluxPostprocessor<2> bx_postprocessor(0, 0);
    MagneticFluxPostprocessor<2> by_postprocessor(0, 1);
    MagneticFluxPostprocessor<2> b_abs_postprocessor(0);
    MatIDPostprocessor<2> mat_id_postprocessor;

    ExportVtu<2> export_vtu(solver.get_triangulation(), solver.get_rhs(), solver.get_solution(), solver.get_fe());
    export_vtu.attach_postprocessor(&mat_id_postprocessor, "MatID");
    export_vtu.attach_postprocessor(&b_abs_postprocessor, "|B| [T]");
    export_vtu.attach_postprocessor(&bx_postprocessor, "Bx [T]");
    export_vtu.attach_postprocessor(&by_postprocessor, "By [T]");
    export_vtu.write("amp");

    return 0;
}