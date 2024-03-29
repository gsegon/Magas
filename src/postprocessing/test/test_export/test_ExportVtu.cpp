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
#include <unordered_map>
#include <string>
#include <any>
#include <fstream>
#include <filesystem>
#include <utility>
#include <variant>

#include "LinearSolver.h"
#include "export/ExportVtu.h"
#include "processors/MatIDPostprocessor.h"
#include "NuCurve.h"
#include "LinearNuCurve.h"
#include "ConstFSource.h"

TEST(ExportVtu, initialize_unit){

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


    std::filesystem::path home = std::getenv("HOME");
    std::filesystem::path test_mesh = "../../../examples/EI_core/EI_core.msh";


    std::unordered_map<int, NuCurve*> nu_map{{200, new LinearNuCurve{nu_core}},       // Core1
                                             {201, new LinearNuCurve{nu_core}},       // Core2
                                             {202, new LinearNuCurve{nu_0}},       // Copper
                                             {203, new LinearNuCurve{nu_0}},       // Copper
                                             {204, new LinearNuCurve{nu_0}},       // Air
                                             {205, new LinearNuCurve{nu_0}},       // Air
    };

    std::unordered_map<int, std::variant<FSource*, std::pair<double, double>>> f_map{ {200, new ConstFSource{0}},        // Core1
                                           {201, new ConstFSource{0}},        // Core2
                                           {202, new ConstFSource{J1}},       // Copper
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

    // Export

    MatIDPostprocessor<2> mat_id_postprocessor;

    ExportVtu<2> export_vtu(solver.get_triangulation(), solver.get_rhs(), solver.get_solution(), solver.get_fe());
    export_vtu.attach_postprocessor(&mat_id_postprocessor, "MatID");

    export_vtu.write("vtu_export_EI_core");

}


