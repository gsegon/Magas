//
// Created by gordan on 3/12/23.
//
#include <gtest/gtest.h>
#include <tuple>
#include <unordered_map>
#include <string>
#include <fstream>
#include <variant>
#include <filesystem>

#include "LinearSolver.h"
#include "NuCurve.h"
#include "LinearNuCurve.h"


using namespace dealii;

TEST(LinearSolver, instantiation){

    LinearSolver<2> solver;

}

TEST(LinearSolver, read_mesh){

    std::filesystem::path home = std::getenv("HOME");
    std::filesystem::path test_mesh = "../../../examples/unit_square/unit_square.msh";


    LinearSolver<2> solver;
    solver.read_mesh(test_mesh);

}

TEST(LinearSolver, setup_system){

    std::filesystem::path home = std::getenv("HOME");
    std::filesystem::path test_mesh = "../../../examples/unit_square/unit_square.msh";


    LinearSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.setup_system();

}

TEST(LinearSolver, assemble_system){

    std::filesystem::path home = std::getenv("HOME");
    std::filesystem::path test_mesh = "../../../examples/unit_square/unit_square.msh";

    std::unordered_map<int, NuCurve*> nu_map{{6, new LinearNuCurve{1}}};
    std::unordered_map<int, std::variant<double, std::pair<double, double>>> f_map{{6, 1}};
    std::unordered_map<int, double> dc_map{{5, 0}};

    LinearSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.set_nu_map(nu_map);
    solver.set_f_map(f_map);
    solver.set_dc_map(dc_map);
    solver.setup_system();
    solver.assemble_system();

}

TEST(LinearSolver, solve_system){

    std::filesystem::path home = std::getenv("HOME");
    std::filesystem::path test_mesh = "../../../examples/unit_square/unit_square.msh";

    std::unordered_map<int, NuCurve*> nu_map{{6, new LinearNuCurve{1}}};
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

}

TEST(LinearSolver, set_nu_map){

    std::unordered_map<int, NuCurve*> nu_map{{6, new LinearNuCurve{1}}};
    std::unordered_map<int, std::variant<double, std::pair<double, double>>> f_map{{6, 1}};
    std::unordered_map<int, double> dc_map{{5, 0}};

    LinearSolver<2> solver;
    solver.set_nu_map(nu_map);
}

TEST(LinearSolver, set_f_map){

    std::unordered_map<int, NuCurve*> nu_map{{6, new LinearNuCurve{1}}};
    std::unordered_map<int, std::variant<double, std::pair<double, double>>> f_map{{6, 1}};
    std::unordered_map<int, double> dc_map{{5, 0}};

    LinearSolver<2> solver;
    solver.set_f_map(f_map);
}

TEST(LinearSolver, setup_system_3){

    std::filesystem::path home = std::getenv("HOME");
    std::filesystem::path test_mesh = "../../../examples/2_conductors/2_conductors_x.msh";


    LinearSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.setup_system();

}

TEST(LinearSolver, 2_conductors){

    constexpr double mu_0 = 1.2566370621219e-6;
    constexpr double nu_0 = 1/mu_0;
    double i_current = 1e3;
    double Jdensity = i_current / (std::pow(0.1,2) * M_PI);


    std::filesystem::path home = std::getenv("HOME");
    std::filesystem::path test_mesh = "../../../examples/2_conductors/2_conductors_x.msh";

    std::unordered_map<int, NuCurve*> nu_map{{1, new LinearNuCurve{nu_0}},       // Conductor 1
                                           {  2, new LinearNuCurve{nu_0}},       // Conductor 2
                                           {  3, new LinearNuCurve{nu_0}}};      // Air

    std::unordered_map<int, std::variant<double, std::pair<double, double>>> f_map{ {1, Jdensity},   // Conductor 1
                                           {2, -Jdensity},  // Conductor 2
                                           {3, 0},          // Air
                                        };

    std::unordered_map<int, double> dc_map{{100, 0}};       // Outerbounds

    LinearSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.set_nu_map(nu_map);
    solver.set_f_map(f_map);
    solver.set_dc_map(dc_map);
    solver.setup_system();
    solver.assemble_system();
    solver.solve();

}

TEST(LinearSolver, EI_core){

    constexpr double mu_0 = 1.2566370621219e-6;
    constexpr double nu_0 = 1/mu_0;
    constexpr double nu_core = 1/(2500*mu_0);
    double J1 = 10*66/8.0645e-05;
    double J2 = -10*66/8.0645e-05;

    std::filesystem::path home = std::getenv("HOME");
    std::filesystem::path test_mesh = "../../../examples/EI_core/EI_core.msh";

    std::unordered_map<int, NuCurve*> nu_map{{200, new LinearNuCurve{nu_core}},       // Core1
                                             {201, new LinearNuCurve{nu_core}},       // Core2
                                             {202, new LinearNuCurve{nu_0}},       // Copper
                                             {203, new LinearNuCurve{nu_0}},       // Copper
                                             {204, new LinearNuCurve{nu_0}},       // Air
                                             {205, new LinearNuCurve{nu_0}},       // Air
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

}

TEST(LinearSolver, motoric_section){

    constexpr double mu_0 = 1.2566370621219e-6;
    constexpr double nu_0 = 1/mu_0;
    constexpr double nu_core = 1/(2500*mu_0);
    double J1 = 10*66/8.0645e-05;
    double J2 = -10*66/8.0645e-05;

    std::filesystem::path home = std::getenv("HOME");
    std::filesystem::path test_mesh = "../../../examples/motoric_section/motoric_section.msh";


    std::unordered_map<int, NuCurve*> nu_map{{1,   new LinearNuCurve{nu_core}},            // Core rotor
                                           {  2,   new LinearNuCurve{nu_core}},            // Core stator

                                           {  506, new LinearNuCurve{nu_0}},             // copper
                                           {  507, new LinearNuCurve{nu_0}},
                                             {508, new LinearNuCurve{nu_0}},
                                             {509, new LinearNuCurve{nu_0}},
                                             {510, new LinearNuCurve{nu_0}},
                                             {511, new LinearNuCurve{nu_0}},

                                             {512, new LinearNuCurve{nu_0}},         // magnets
                                           {513, new LinearNuCurve{nu_0}},
                                             {514, new LinearNuCurve{nu_0}},
                                             {515, new LinearNuCurve{nu_0}},

                                             {516, new LinearNuCurve{nu_0}},          // Air
                                           {517, new LinearNuCurve{nu_0}},
    };

    std::unordered_map<int, std::variant<double, std::pair<double, double>>> f_map{ {1, 0},       // Core rotor
                                                                                    {2, 0},       // Core stator

                                                                                    {506, 0},         // Coils
                                                                                    {507, 0},
                                                                                    {508, J1},
                                                                                    {509, J1},
                                                                                    {510, 0},
                                                                                    {511, 0},

                                                                                    {512, 0},         // magnets
                                                                                    {513, 0},
                                                                                    {514, 0},
                                                                                    {515, 0},
                                                                                    {516, 0},
                                                                                    {517, 0},
    };

    std::unordered_map<int, double> dc_map{{505, 0.0},
//                                           {518, 0.0},
//                                           {521, 0.0},
//                                           {519, 0.0},
//                                           {520, 0.0}
    } ;

    // Solve
    LinearSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.set_nu_map(nu_map);
    solver.set_f_map(f_map);
    solver.set_dc_map(dc_map);
    solver.setup_system();
    solver.assemble_system();
    solver.solve();

}

TEST(LinearSolver, Magnet){

    constexpr double mu_0 = 1.2566370621219e-6;
    constexpr double nu_0 = 1/mu_0;

    std::filesystem::path home = std::getenv("HOME");
    std::filesystem::path test_mesh = "../../../examples/magnet/BlockMagnet.msh";


    std::unordered_map<int, NuCurve*> nu_map{{1, new LinearNuCurve{nu_0}},       // Air
                                           {  2, new LinearNuCurve{nu_0}},       // Magnet

    };

    std::pair<double, double> Hc{0, 10e3};
    std::unordered_map<int, std::variant<double, std::pair<double, double>>> f_map{ {1, 0},
                                           {2, Hc},        // Magnet
                                           };

    std::unordered_map<int, double> dc_map{{505, 0}} ;      // Outer boundary

    // Solve
    LinearSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.set_nu_map(nu_map);
    solver.set_f_map(f_map);
    solver.set_dc_map(dc_map);
    solver.setup_system();
    solver.assemble_system();
    solver.solve();

}


