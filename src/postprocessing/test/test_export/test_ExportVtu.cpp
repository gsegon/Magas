//
// Created by gordan on 3/12/23.
//
#include <gtest/gtest.h>
#include <unordered_map>
#include <string>
#include <any>
#include <fstream>

#include "LinearSolver.h"
#include "export/ExportVtu.h"
#include "processors/MatIDPostprocessor.h"



TEST(ExportVtu, initialize_unit){

    std::string test_mesh = "~/Programs/solver/examples/unit_square/unit_square.msh";
    std::unordered_map<int, double> nu_map{{6, 1}};
    std::unordered_map<int, std::variant<double, std::pair<double, double>>> f_map{{6, 1}};
    std::unordered_map<int, double> dc_map{{5, 0}};

    LinearSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.set_nu_map(nu_map);
    solver.set_f_map(f_map);
    solver.set_dc_map(dc_map);
    solver.setup_system();
    solver.assemble_system();
    solver.solve();

    // Export
    ExportVtu<2> export_vtu3(solver.get_triangulation(), solver.get_rhs(), solver.get_solution(), solver.get_fe());
    export_vtu3.write("vtu_export_unit_square_core");

}


TEST(ExportVtu, initialize_EI_core){

    constexpr double mu_0 = 1.2566370621219e-6;
    constexpr double nu_0 = 1/mu_0;
    constexpr double nu_core = 1/(2500*mu_0);
    double J1 = 10*66/8.0645e-05;
    double J2 = -10*66/8.0645e-05;


    std::string test_mesh = "~/Programs/solver/examples/EI_core/EI_core.msh";
    std::unordered_map<int, double> nu_map{{200, nu_core},      // Core1
                                           {201, nu_core},      // Core2
                                           {202, nu_0},         // Copper
                                           {203, nu_0},         // Copper
                                           {204, nu_0},         // Air
                                           {205, nu_0},         // Air
    };

    std::unordered_map<int, std::variant<double, std::pair<double, double>>> f_map{ {200, 0},        // Core1
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
    solver.set_nu_map(nu_map);
    solver.set_f_map(f_map);
    solver.set_dc_map(dc_map);
    solver.setup_system();
    solver.assemble_system();
    solver.solve();

    // Export

    MatIDPostprocessor<2> mat_id_postprocessor;

    ExportVtu<2> export_vtu(solver.get_triangulation(), solver.get_rhs(), solver.get_solution(), solver.get_fe());
    export_vtu.attach_postprocessor(&mat_id_postprocessor, "MatID");

    export_vtu.write("vtu_export_EI_core");

}


