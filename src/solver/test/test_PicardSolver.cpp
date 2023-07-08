//
// Created by gordan on 3/12/23.
//
#include <gtest/gtest.h>
#include <tuple>
#include <unordered_map>
#include <string>
#include <any>
#include <fstream>
#include <filesystem>

#include "PicardSolver.h"

TEST(PicardSolver, instantiation){

    PicardSolver<2> solver;

}

TEST(PicardSolver, read_mesh){

    std::filesystem::path home = std::getenv("HOME");
    std::filesystem::path test_mesh = "../../../examples/EI_core/EI_core.msh";

    PicardSolver<2> solver;
    solver.read_mesh(test_mesh);


}

TEST(PicardSolver, setup_cell_nu_history){

    std::filesystem::path home = std::getenv("HOME");
    std::filesystem::path test_mesh = "../../../examples/EI_core/EI_core.msh";

    PicardSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.setup_cell_nu_history();

}

TEST(PicardSolver, setup_system){

    std::filesystem::path home = std::getenv("HOME");
    std::filesystem::path test_mesh = "../../../examples/EI_core/EI_core.msh";

    PicardSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.setup_cell_nu_history();
    solver.setup_system();

}

TEST(PicardSolver, reinit_system){

    std::filesystem::path home = std::getenv("HOME");
    std::filesystem::path test_mesh = "../../../examples/unit_square/unit_square.msh";

    std::unordered_map<int, std::any> nu_map{{6, "Nonlinear"}};
    std::unordered_map<int, double> f_map{{6, 1.0}};
    std::unordered_map<int, double> dc_map{{5, 0.0}};

    PicardSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.setup_cell_nu_history();
    solver.setup_system();
    solver.reinit_system();

}

TEST(PicardSolver, assemble_system){

    std::filesystem::path home = std::getenv("HOME");
    std::filesystem::path test_mesh = "../../../examples/unit_square/unit_square.msh";

    std::unordered_map<int, std::any> nu_map{{6, "Nonlinear"}};
    std::unordered_map<int, double> f_map{{6, 1.0}};
    std::unordered_map<int, double> dc_map{{5, 0.0}};

    PicardSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.setup_cell_nu_history();
    solver.set_nu_map(nu_map);
    solver.set_f_map(f_map);
    solver.set_dc_map(dc_map);
    solver.setup_system();

    solver.reinit_system();
    solver.assemble_system();

}

TEST(PicardSolver, solve_system){

    std::filesystem::path home = std::getenv("HOME");
    std::filesystem::path test_mesh = "../../../examples/unit_square/unit_square.msh";

    std::unordered_map<int, std::any> nu_map{{6, "Nonlinear"}};
    std::unordered_map<int, double> f_map{{6, 1.0}};
    std::unordered_map<int, double> dc_map{{5, 0.0}};

    PicardSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.setup_cell_nu_history();
    solver.set_nu_map(nu_map);
    solver.set_f_map(f_map);
    solver.set_dc_map(dc_map);
    solver.setup_system();
    solver.initialize_cell_nu_history(1);

    solver.reinit_system();
    solver.assemble_system();
    solver.solve();

}



TEST(PicardSolver, solve_nonlinear){

    std::filesystem::path home = std::getenv("HOME");
    std::filesystem::path test_mesh = "../../../examples/unit_square/unit_square.msh";

    std::unordered_map<int, std::any> nu_map{{6, "Nonlinear"}};
    std::unordered_map<int, double> f_map{{6, 1.0}};
    std::unordered_map<int, double> dc_map{{5, 0.0}};

    PicardSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.setup_cell_nu_history();
    solver.set_nu_map(nu_map);
    solver.set_f_map(f_map);
    solver.set_dc_map(dc_map);
    solver.setup_system();
    solver.initialize_cell_nu_history(1);
    solver.solve_nonlinear(3);

}


TEST(PicardSolver, EI_core){

    constexpr double mu_0 = 1.2566370621219e-6;
    constexpr double nu_0 = 1/mu_0;
    constexpr double nu_core = 1/(2500*mu_0);
    double J1 = 30*66/8.0645e-05;
    double J2 = -30*66/8.0645e-05;


    std::filesystem::path home = std::getenv("HOME");
    std::filesystem::path test_mesh = "../../../examples/EI_core/EI_core.msh";

    std::unordered_map<int, std::any> nu_map{{200, "Nonlinear"},       // Core1
                                             {201, "Nonlinear"},       // Core2
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

    PicardSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.setup_cell_nu_history();
    solver.set_nu_map(nu_map);
    solver.set_f_map(f_map);
    solver.set_dc_map(dc_map);
    solver.setup_system();
    solver.reinit_system();
    solver.initialize_cell_nu_history(nu_core);
    solver.solve_nonlinear(3);

}
