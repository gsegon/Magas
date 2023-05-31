//
// Created by gordan on 3/12/23.
//
#include <gtest/gtest.h>
#include <tuple>
#include <unordered_map>
#include <string>
#include <fstream>
#include <variant>

#include "NewtonSolver.h"
#include "MagneticFluxPostprocessor.h"
#include "MagneticEnergyPostprocessor.h"
#include "MagneticEnergyDensityPostprocessor.h"
#include "ExportVtu.h"
#include "MatIDPostprocessor.h"
#include <deal.II/numerics/data_out.h>

using namespace dealii;

TEST(NewtonSolver, instantiation){

    NewtonSolver<2> solver;

}

TEST(NewtonSolver, read_mesh){
    std::string test_mesh = "/home/gordan/Programs/solver/test/test_data/test_unit_square/unit_square.msh";

    NewtonSolver<2> solver;
    solver.read_mesh(test_mesh);


}

TEST(NewtonSolver, setup_system){

    std::string test_mesh = "/home/gordan/Programs/solver/test/test_data/test_unit_square/unit_square.msh";

    NewtonSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.setup_system(true);

}

TEST(NewtonSolver, set_maps){

    std::string test_mesh = "/home/gordan/Programs/solver/test/test_data/test_unit_square/unit_square.msh";
    std::unordered_map<int, std::any> nu_map{{6, 1}};
    std::unordered_map<int, std::variant<double, std::pair<double, double>>> f_map{{6, 1}};
    std::unordered_map<int, double> dc_map{{5, 0}};

    NewtonSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.setup_system(true);
    solver.set_nu_map(nu_map);
    solver.set_f_map(f_map);
    solver.set_dc_map(dc_map);


}
TEST(NewtonSolver, assemble){

    std::string test_mesh = "/home/gordan/Programs/solver/test/test_data/test_unit_square/unit_square.msh";
    std::unordered_map<int, std::any> nu_map{{6, 1}};
    std::unordered_map<int, std::variant<double, std::pair<double, double>>> f_map{{6, 1}};
    std::unordered_map<int, double> dc_map{{5, 0}};

    NewtonSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.setup_system(true);
    solver.set_nu_map(nu_map);
    solver.set_f_map(f_map);
    solver.set_dc_map(dc_map);
    solver.assemble_system();

}

TEST(NewtonSolver, solve){

    std::string test_mesh = "/home/gordan/Programs/solver/test/test_data/test_unit_square/unit_square.msh";
    std::unordered_map<int, std::any> nu_map{{6, 1.0}};
    std::unordered_map<int, std::variant<double, std::pair<double, double>>> f_map{{6, 1}};
    std::unordered_map<int, double> dc_map{{5, 0}};

    NewtonSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.setup_system(true);
    solver.set_nu_map(nu_map);
    solver.set_f_map(f_map);
    solver.set_dc_map(dc_map);

    ExportVtu<2> export_vtu(solver.get_triangulation(), solver.get_rhs(), solver.get_current_solution(), solver.get_fe());

    int i = 0;
    double alpha = 0.0;
    while(i++ < 50){
        solver.assemble_system();
        if (i < 5){
            alpha = 0.1;
        }
        else if (i < 15){
            alpha = 0.3;
        }
        else{
            alpha = 0.5;
        }

        solver.solve(alpha);
        std::cout << "\tResidual(" << i<< "): " << solver.compute_residual() << std::endl;
        export_vtu.write("vtu_export_newton_unit_square-" + std::to_string(i));
    }


}

TEST(NewtonSolver, EI_core){

    constexpr double mu_0 = 1.2566370621219e-6;
    constexpr double nu_0 = 1/mu_0;
    constexpr double nu_core = 1/(2500*mu_0);
    double J1 = 1*66/8.0645e-05;
    double J2 = -1*66/8.0645e-05;


    std::string test_mesh = "/home/gordan/Programs/solver/test/test_data/test_EI_core/EI_core.msh";
    std::unordered_map<int, std::any> nu_map{{200, "Nonlinear"},       // Core1
                                             {201, "Nonlinear"},       // Core2
                                             {202, nu_0},       // Copper
                                             {203, nu_0},       // Copper
                                             {204, nu_0},       // Air
                                             {205, nu_0},       // Air
    };

    std::unordered_map<int, std::variant<double, std::pair<double, double>>> f_map{ {200, 0.0},        // Core1
                                           {201, 0.0},        // Core2
                                           {202, J1},       // Copper
                                           {203, J2},       // Copper
                                           {204, 0.0},        // Air
                                           {205, 0.0},        // Air
    };

    std::unordered_map<int, double> dc_map{{44, 0}} ;      // Air

    NewtonSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.setup_system(true);
    solver.set_nu_map(nu_map);
    solver.set_f_map(f_map);
    solver.set_dc_map(dc_map);

    ExportVtu<2> export_vtu(solver.get_triangulation(), solver.get_rhs(), solver.get_current_solution(), solver.get_fe());

    // Export post-processors
    MagneticFluxPostprocessor<2> b_abs_postprocessor_q0(0);
    MagneticFluxPostprocessor<2> b_abs_postprocessor_q1(1);
    MagneticFluxPostprocessor<2> b_abs_postprocessor_q2(2);
    MagneticFluxPostprocessor<2> b_abs_postprocessor_q3(3);


    export_vtu.attach_postprocessor(&b_abs_postprocessor_q0, "Babs_q0 [T]");
    export_vtu.attach_postprocessor(&b_abs_postprocessor_q1, "Babs_q1 [T]");
    export_vtu.attach_postprocessor(&b_abs_postprocessor_q2, "Babs_q2 [T]");
    export_vtu.attach_postprocessor(&b_abs_postprocessor_q3, "Babs_q3 [T]");

    int i = 0;
    double alpha = 0.0;
    while(i++ < 50){

        if (i < 10)
            alpha = 0.1;
        else if (i < 15)
            alpha = 0.3;
        else
            alpha = 0.5;

        solver.assemble_system();
        solver.solve(alpha);
        double res = solver.compute_residual();
        std::cout << "\tResidual(" << i<< "): " << res << std::endl;
        export_vtu.write("vtu_export_newton_EI_core-" + std::to_string(i));
        if (res < 1e-6){
            std::cout << "Converged!";
            break;
        }

    }

}
//
//TEST(NewtonSolver, Magnet){
//
//    constexpr double mu_0 = 1.2566370621219e-6;
//    constexpr double nu_0 = 1/mu_0;
//
//
//    std::string test_mesh = "/home/gordan/Programs/solver/test/test_data/test_magnet/BlockMagnet.msh";
//    std::unordered_map<int, double> nu_map{{1, nu_0},       // Air
//                                           {2, nu_0},       // Magnet
//
//    };
//
//    std::pair<double, double> Hc{0, 10e3};
//    std::unordered_map<int, std::variant<double, std::pair<double, double>>> f_map{ {1, 0},
//                                           {2, Hc},        // Magnet
//                                           };
//
//    std::unordered_map<int, double> dc_map{{505, 0}} ;      // Outer boundary
//
//    // Solve
//    NewtonSolver<2> solver;
//    solver.read_mesh(test_mesh);
//    solver.setup_system();
//    solver.set_nu_map(nu_map);
//    solver.set_f_map(f_map);
//    solver.set_dc_map(dc_map);
//    solver.assemble_system();
//    solver.solve();
//
//    // Export
//    MagneticFluxPostprocessor<2> bx_postprocessor_q0(0, 0);
//    MagneticFluxPostprocessor<2> by_postprocessor_q0(0, 1);
//    MagneticFluxPostprocessor<2> bx_postprocessor_q1(1, 0);
//    MagneticFluxPostprocessor<2> by_postprocessor_q1(1, 1);
//    MagneticFluxPostprocessor<2> bx_postprocessor_q2(2, 0);
//    MagneticFluxPostprocessor<2> by_postprocessor_q2(2, 1);
//    MagneticFluxPostprocessor<2> bx_postprocessor_q3(3, 0);
//    MagneticFluxPostprocessor<2> by_postprocessor_q3(3, 1);
//
//    MagneticEnergyPostprocessor<2> energy_cell(nu_map);
//    MagneticEnergyDensityPostprocessor<2> energy_density(nu_map);
//
//    MatIDPostprocessor<2> mat_id_postprocessor;
//
//    ExportVtu<2> export_vtu(solver.get_triangulation(), solver.get_rhs(), solver.get_solution(), solver.get_fe());
//    export_vtu.attach_postprocessor(&mat_id_postprocessor, "MatID");
//
//    export_vtu.attach_postprocessor(&bx_postprocessor_q0, "Bx_q0 [T]");
//    export_vtu.attach_postprocessor(&by_postprocessor_q0, "By_q0 [T]");
//    export_vtu.attach_postprocessor(&bx_postprocessor_q1, "Bx_q1 [T]");
//    export_vtu.attach_postprocessor(&by_postprocessor_q1, "By_q1 [T]");
//    export_vtu.attach_postprocessor(&bx_postprocessor_q2, "Bx_q2 [T]");
//    export_vtu.attach_postprocessor(&by_postprocessor_q2, "By_q2 [T]");
//    export_vtu.attach_postprocessor(&bx_postprocessor_q3, "Bx_q3 [T]");
//    export_vtu.attach_postprocessor(&by_postprocessor_q3, "By_q3 [T]");
//
//    export_vtu.attach_postprocessor(&energy_cell, "E [J/m]");
//    export_vtu.attach_postprocessor(&energy_density, "E [J/m^3]");
//
//    export_vtu.write("vtu_export_block_magnet");
//
//}
//
