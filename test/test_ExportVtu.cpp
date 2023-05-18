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
#include "MagneticFluxPostprocessor.h"
#include "MatIDPostprocessor.h"
#include "MagneticEnergyPostprocessor.h"
#include "MagneticEnergyDensityPostprocessor.h"



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
    export_vtu3.write("vtu_export_unit_square_core");

}


TEST(ExportVtu, initialize_EI_core){

    constexpr double mu_0 = 1.2566370621219e-6;
    constexpr double nu_0 = 1/mu_0;
    constexpr double nu_core = 1/(2500*mu_0);
    double J1 = 10*66/8.0645e-05;
    double J2 = -10*66/8.0645e-05;


    std::string test_mesh = "/home/gordan/Programs/solver/test/test_data/test_EI_core/EI_core.msh";
    std::unordered_map<int, double> nu_map{{200, nu_core},      // Core1
                                           {201, nu_core},      // Core2
                                           {202, nu_0},         // Copper
                                           {203, nu_0},         // Copper
                                           {204, nu_0},         // Air
                                           {205, nu_0},         // Air
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
    MagneticFluxPostprocessor<2> bx_postprocessor_q0(0, 0);
    MagneticFluxPostprocessor<2> by_postprocessor_q0(0, 1);
    MagneticFluxPostprocessor<2> bx_postprocessor_q1(1, 0);
    MagneticFluxPostprocessor<2> by_postprocessor_q1(1, 1);
    MagneticFluxPostprocessor<2> bx_postprocessor_q2(2, 0);
    MagneticFluxPostprocessor<2> by_postprocessor_q2(2, 1);
    MagneticFluxPostprocessor<2> bx_postprocessor_q3(3, 0);
    MagneticFluxPostprocessor<2> by_postprocessor_q3(3, 1);

    MagneticEnergyPostprocessor<2> energy_cell(nu_map);
    MagneticEnergyDensityPostprocessor<2> energy_density(nu_map);

    MatIDPostprocessor<2> mat_id_postprocessor;

    ExportVtu<2> export_vtu(solver.get_triangulation(), solver.get_rhs(), solver.get_solution(), solver.get_fe());
    export_vtu.attach_postprocessor(&mat_id_postprocessor, "MatID");

    export_vtu.attach_postprocessor(&bx_postprocessor_q0, "Bx_q0 [T]");
    export_vtu.attach_postprocessor(&by_postprocessor_q0, "By_q0 [T]");
    export_vtu.attach_postprocessor(&bx_postprocessor_q1, "Bx_q1 [T]");
    export_vtu.attach_postprocessor(&by_postprocessor_q1, "By_q1 [T]");
    export_vtu.attach_postprocessor(&bx_postprocessor_q2, "Bx_q2 [T]");
    export_vtu.attach_postprocessor(&by_postprocessor_q2, "By_q2 [T]");
    export_vtu.attach_postprocessor(&bx_postprocessor_q3, "Bx_q3 [T]");
    export_vtu.attach_postprocessor(&by_postprocessor_q3, "By_q3 [T]");

    export_vtu.attach_postprocessor(&energy_cell, "E [J/m]");
    export_vtu.attach_postprocessor(&energy_density, "E [J/m^3]");

    export_vtu.write("vtu_export_EI_core");

}


