//
// Created by gordan on 5/12/23.
//

#include <iostream>
#include <fstream>
#include "LinearSolver.h"
#include "nlohmann/json.hpp"

#include "ExportVtu.h"
#include "MagneticFluxPostprocessor.h"
#include "MatIDPostprocessor.h"
#include "MagneticEnergyDensityPostprocessor.h"
#include "MagneticEnergyPostprocessor.h"

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

    std::unordered_map<int, std::variant<double, std::pair<double, double>>> f_map;
    std::unordered_map<int, double> nu_map;
    std::unordered_map<int, double> dc_map;

    for (auto& mesh_el_data : mesh_id_data.items()){
        int mat_id{std::stoi(mesh_el_data.key())};
        if (mesh_el_data.value().contains("material")){
            nu_map.insert({mat_id, material_data.at(mesh_el_data.value().at("material")).at("nu")});
        }
        if (mesh_el_data.value().contains("source")){

            auto source_d = source_data.at(mesh_el_data.value().at("source"));

            // If mesh data contains angle, the source is vector Hc. Magnitude in Material data and direction from mesh data.
            if (mesh_el_data.value().contains("angle")){
                double Hc = source_d;
                double angle = mesh_el_data.value().at("angle");
                double Hc_x = Hc*std::cos(angle);
                double Hc_y = Hc*std::sin(angle);

                std::pair<double, double> Hc_vec{Hc_x, Hc_y};
                f_map.insert({mat_id, Hc_vec});
            }

            // Otherwise, the source is J (current density)
            else {
                double J = source_d;
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
    export_vtu.write("aml");

    return 0;
}