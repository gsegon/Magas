// Magas - Magnetostatic Analysis Suite
// Copyright (C) 2023  Gordan Segon
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
// USA

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
#include "ConstFSource.h"


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
    std::unordered_map<int, std::variant<FSource*, std::pair<double, double>>> f_map{{6, new ConstFSource{1}}};
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
    std::unordered_map<int, std::variant<FSource*, std::pair<double, double>>> f_map{{6, new ConstFSource{1}}};
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
    std::unordered_map<int, std::variant<FSource*, std::pair<double, double>>> f_map{{6, new ConstFSource{1}}};
    std::unordered_map<int, double> dc_map{{5, 0}};

    LinearSolver<2> solver;
    solver.set_nu_map(nu_map);
}

TEST(LinearSolver, set_f_map){

    std::unordered_map<int, NuCurve*> nu_map{{6, new LinearNuCurve{1}}};
    std::unordered_map<int, std::variant<FSource*, std::pair<double, double>>> f_map{{6, new ConstFSource{1}}};
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

    std::unordered_map<int, std::variant<FSource*, std::pair<double, double>>> f_map{ {1, new ConstFSource{Jdensity}},   // Conductor 1
                                           {2,  new ConstFSource{-Jdensity}},  // Conductor 2
                                           {3,  new ConstFSource{0}},          // Air
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

    std::unordered_map<int, std::variant<FSource*, std::pair<double, double>>> f_map{ {200,  new ConstFSource{0}},        // Core1
                                           {201,  new ConstFSource{0}},        // Core2
                                           {202,  new ConstFSource{J1}},       // Copper
                                           {203, new ConstFSource{J2}},       // Copper
                                           {204, new ConstFSource{0}},        // Air
                                           {205, new ConstFSource{0}},        // Air
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
    std::filesystem::path test_mesh = "../../../examples/ipm_2_section/ipm_2_section.msh";


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

    std::unordered_map<int, std::variant<FSource*, std::pair<double, double>>> f_map{ {1, new ConstFSource{0}},       // Core rotor
                                                                                    {2, new ConstFSource{0}},       // Core stator

                                                                                    {506, new ConstFSource{0}},         // Coils
                                                                                    {507, new ConstFSource{0}},
                                                                                    {508, new ConstFSource{J1}},
                                                                                    {509, new ConstFSource{J2}},
                                                                                    {510, new ConstFSource{0}},
                                                                                    {511, new ConstFSource{0}},

                                                                                    {512, new ConstFSource{0}},         // magnets
                                                                                    {513, new ConstFSource{0}},
                                                                                    {514, new ConstFSource{0}},
                                                                                    {515, new ConstFSource{0}},
                                                                                    {516, new ConstFSource{0}},
                                                                                    {517, new ConstFSource{0}},
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
    std::unordered_map<int, std::variant<FSource*, std::pair<double, double>>> f_map{ {1, new ConstFSource{0}},
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


