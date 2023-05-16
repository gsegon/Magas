//
// Created by gordan on 3/12/23.
//
#include <gtest/gtest.h>
#include <unordered_map>
#include <string>
#include <any>
#include <fstream>

#include "LinearSolver.h"
#include "ExportVtu.h"



TEST(ExportVtu, initialize_unit){

    std::string test_mesh = "/home/gordan/Programs/solver/test/test_data/test_unit_square/unit_square.msh";
    std::unordered_map<int, double> nu_map{{6, 1}};
    std::unordered_map<int, double> f_map{{6, 1}};
    std::unordered_map<int, double> dc_map{{5, 0}};

    LinearSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.setup_system();
    solver.set_nu_map(nu_map);
    solver.set_f_map(f_map);
    solver.set_dc_map(dc_map);
    solver.assemble_system();
    solver.solve();

    // Export
    ExportVtu<2> export_vtu3(solver.get_triangulation(), solver.get_rhs(), solver.get_solution(), solver.get_fe());
    export_vtu3.write("unit_square_core_vtu_export");

}


TEST(ExportVtu, initialize_EI_core){

    constexpr double mu_0 = 1.2566370621219e-6;
    constexpr double nu_0 = 1/mu_0;
    constexpr double nu_core = 1/(2500*mu_0);
    double J1 = 10*66/8.0645e-05;
    double J2 = -10*66/8.0645e-05;


    std::string test_mesh = "/home/gordan/Programs/solver/test/test_data/test_EI_core/EI_core.msh";
    std::unordered_map<int, double> nu_map{{200, nu_core},       // Core1
                                           {201, nu_core},       // Core2
                                           {202, nu_0},       // Copper
                                           {203, nu_0},       // Copper
                                           {204, nu_0},       // Air
                                           {205, nu_0},       // Air
    };

    std::unordered_map<int, double> f_map{ {200, 0},        // Core1
                                           {201, 0},        // Core2
                                           {202, J1},       // Copper
                                           {203, J2},       // Copper
                                           {204, 0},        // Air
                                           {205, 0},        // Air
    };

    std::unordered_map<int, double> dc_map{{44, 0}} ;      // Air

    // Solve
    LinearSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.setup_system();
    solver.set_nu_map(nu_map);
    solver.set_f_map(f_map);
    solver.set_dc_map(dc_map);
    solver.assemble_system();
    solver.solve();

    // Export
    ExportVtu<2> export_vtu3(solver.get_triangulation(), solver.get_rhs(), solver.get_solution(), solver.get_fe());
    export_vtu3.write("EI_core_vtu_export");


}


