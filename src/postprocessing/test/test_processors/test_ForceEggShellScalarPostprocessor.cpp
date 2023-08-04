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
#include <any>
#include <fstream>
#include <filesystem>

#include "LinearSolver.h"
#include "NewtonSolver.h"
#include "processors/ForceEggShellScalarPostprocessor.h"
#include "processors/ExpressionScalarPostprocessor.h"
#include "NuCurve.h"
#include "LinearNuCurve.h"
#include "ConstFSource.h"

TEST(TestForceEggShellScalarPostprocessor, force_benchmark_kelvin_1){

    std::filesystem::path home = std::getenv("HOME");
    std::filesystem::path test_mesh = "../../../examples/force_benchmark_kelvin_1/force_benchmark_kelvin_1.msh";

    std::unordered_map<int, NuCurve*> nu_map{{3, new LinearNuCurve{795774.715025}},
                                             {4, new LinearNuCurve{795774.715025}},
                                             {5, new LinearNuCurve{795774.715025}},
                                             {6, new LinearNuCurve{795774.715025}},
                                             {8, new LinearNuCurve{795774.715025}}};

    std::unordered_map<int, std::variant<FSource*, std::pair<double, double>>> f_map{{3, new ConstFSource{25/(std::pow(0.25e-2, 2)*M_PI)}},
                                                                                   {4, new ConstFSource{0}},
                                                                                   {5, new ConstFSource{0}},
                                                                                   {6, new ConstFSource{0}},
                                                                                   {8, std::pair<double, double>{-1591549.43091895, 0.0}}
    };
    std::unordered_map<int, double> dc_map{{7, 0}};
    std::unordered_map<std::string, std::vector<unsigned int>> per_map{{"periodic-circle", {1, 2}}};

    LinearSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.set_nu_map(nu_map);
    solver.set_f_map(f_map);
    solver.set_dc_map(dc_map);
    solver.set_per_map(per_map);
    solver.setup_system();
    solver.assemble_system();
    solver.solve();

    ForceEggShellScalarPostprocessor<2> force_x_eggshell_postp{3, 4, "x"};
    ForceEggShellScalarPostprocessor<2> force_y_eggshell_postp{3, 4, "y"};
    ExpressionScalarPostprocessor<2> force_x_lorentz_postp{"if((mat_id == 3), -4e-2*(By_q1*J_q1*JxW_q1+By_q2*J_q2*JxW_q2+By_q3*J_q3*JxW_q3+By_q4*J_q4*JxW_q4), 0)",
                                                           nu_map,
                                                           f_map};
    ExpressionScalarPostprocessor<2> force_y_lorentz_postp{"if((mat_id == 3), 4e-2*(Bx_q1*J_q1*JxW_q1+Bx_q2*J_q2*JxW_q2+Bx_q3*J_q3*JxW_q3+Bx_q4*J_q4*JxW_q4), 0)",
                                                           nu_map,
                                                           f_map};

    double f_x;
    force_x_eggshell_postp.process(solver.get_triangulation(), solver.get_solution(), solver.get_fe(), f_x);

    double f_y;
    force_y_eggshell_postp.process(solver.get_triangulation(), solver.get_solution(), solver.get_fe(), f_y);

    double f_x_expression;
    force_x_lorentz_postp.process(solver.get_triangulation(), solver.get_solution(), solver.get_fe(), f_x_expression);

    double f_y_expression;
    force_y_lorentz_postp.process(solver.get_triangulation(), solver.get_solution(), solver.get_fe(), f_y_expression);

    std::cout << "f_x: " << 4e-2*f_x << std::endl;
    std::cout << "f_x_expression: " << f_x_expression << std::endl;
    std::cout << "f_y: " << 4e-2*f_y << std::endl;
    std::cout << "f_y_expression: " << f_y_expression << std::endl;

    ASSERT_NEAR(f_y*4e-2, -1.0, 1e-3);
    ASSERT_NEAR(f_x*4e-2, 0, 1e-3);

    ASSERT_NEAR(f_y*4e-2, f_y_expression, 1e-3);
    ASSERT_NEAR(f_x*4e-2, f_x_expression, 1e-3);

}

TEST(TestForceEggShellScalarPostprocessor, force_benchmark_kelvin_2){

    std::filesystem::path home = std::getenv("HOME");
    std::filesystem::path test_mesh = "../../../examples/force_benchmark_kelvin_1/force_benchmark_kelvin_1.msh";

    std::unordered_map<int, NuCurve*> nu_map{{3, new LinearNuCurve{795774.715025}},
                                             {4, new LinearNuCurve{795774.715025}},
                                             {5, new LinearNuCurve{795774.715025}},
                                             {6, new LinearNuCurve{795774.715025}},
                                             {8, new LinearNuCurve{795774.715025}}};

    std::unordered_map<int, std::variant<FSource*, std::pair<double, double>>> f_map{{3, new ConstFSource{25/(std::pow(0.25e-2, 2)*M_PI)}},
                                                                                   {4, new ConstFSource{0}},
                                                                                   {5, new ConstFSource{0}},
                                                                                   {6, new ConstFSource{0}},
                                                                                   {8, std::pair<double, double>{0.0, -1591549.43091895}}
    };
    std::unordered_map<int, double> dc_map{{7, 0}};
    std::unordered_map<std::string, std::vector<unsigned int>> per_map{{"periodic-circle", {1, 2}}};

    LinearSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.set_nu_map(nu_map);
    solver.set_f_map(f_map);
    solver.set_dc_map(dc_map);
    solver.set_per_map(per_map);
    solver.setup_system();
    solver.assemble_system();
    solver.solve();

    ForceEggShellScalarPostprocessor<2> force_x_eggshell_postp{3, 4, "x"};
    ForceEggShellScalarPostprocessor<2> force_y_eggshell_postp{3, 4, "y"};
    ExpressionScalarPostprocessor<2> force_x_lorentz_postp{"if((mat_id == 3), -4e-2*(By_q1*J_q1*JxW_q1+By_q2*J_q2*JxW_q2+By_q3*J_q3*JxW_q3+By_q4*J_q4*JxW_q4), 0)",
                                                           nu_map,
                                                           f_map};
    ExpressionScalarPostprocessor<2> force_y_lorentz_postp{"if((mat_id == 3), 4e-2*(Bx_q1*J_q1*JxW_q1+Bx_q2*J_q2*JxW_q2+Bx_q3*J_q3*JxW_q3+Bx_q4*J_q4*JxW_q4), 0)",
                                                           nu_map,
                                                           f_map};

    double f_x;
    force_x_eggshell_postp.process(solver.get_triangulation(), solver.get_solution(), solver.get_fe(), f_x);

    double f_y;
    force_y_eggshell_postp.process(solver.get_triangulation(), solver.get_solution(), solver.get_fe(), f_y);

    double f_x_expression;
    force_x_lorentz_postp.process(solver.get_triangulation(), solver.get_solution(), solver.get_fe(), f_x_expression);

    double f_y_expression;
    force_y_lorentz_postp.process(solver.get_triangulation(), solver.get_solution(), solver.get_fe(), f_y_expression);

    std::cout << "f_x: " << 4e-2*f_x << std::endl;
    std::cout << "f_x_expression: " << f_x_expression << std::endl;
    std::cout << "f_y: " << 4e-2*f_y << std::endl;
    std::cout << "f_y_expression: " << f_y_expression << std::endl;

    ASSERT_NEAR(f_y*4e-2, 0, 1e-3);
    ASSERT_NEAR(f_x*4e-2, 1, 1e-3);

    ASSERT_NEAR(f_y*4e-2, f_y_expression, 1e-3);
    ASSERT_NEAR(f_x*4e-2, f_x_expression, 1e-3);

}


