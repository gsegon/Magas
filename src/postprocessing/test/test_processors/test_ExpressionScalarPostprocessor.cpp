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
#include "NuCurve.h"
#include "LinearNuCurve.h"
#include "processors/ExpressionScalarPostprocessor.h"
#include "ConstFSource.h"
#include "FSource.h"

TEST(ExpressionScalarPostprocessor, unit_square){

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

    ExpressionScalarPostprocessor<2> expression_postp_count_elem{"1", nu_map};
    ExpressionScalarPostprocessor<2> expression_postp_energy{"(Bx_q1^2+By_q1^2)/2.0 * nu_q1 * JxW_q1 + "
                                                             "(Bx_q2^2+By_q2^2)/2.0 * nu_q2 * JxW_q2 + "
                                                             "(Bx_q3^2+By_q3^2)/2.0 * nu_q3 * JxW_q3 + "
                                                             "(Bx_q4^2+By_q4^2)/2.0 * nu_q4 * JxW_q4", nu_map};

    double result_count_elem;
    expression_postp_count_elem.process(solver.get_triangulation(), solver.get_solution(), solver.get_fe(), result_count_elem);

    double result_energy =0;
    expression_postp_energy.process(solver.get_triangulation(), solver.get_solution(), solver.get_fe(), result_energy);

    std::cout << "result_count_elem: " << result_count_elem << std::endl;
    std::cout << "result_energy: " << result_energy << std::endl;

}

TEST(ExpressionScalarPostprocessorNonlinear, unit_square){

    std::filesystem::path home = std::getenv("HOME");
    std::filesystem::path test_mesh = "../../../examples/unit_square/unit_square.msh";


    std::unordered_map<int, NuCurve*> nu_map{{6, new LinearNuCurve{1}}};
    std::unordered_map<int, std::variant<FSource*, std::pair<double, double>>> f_map{{6, new ConstFSource{1}}};
    std::unordered_map<int, double> dc_map{{5, 0}};

    NewtonSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.set_nu_map(nu_map);
    solver.set_f_map(f_map);
    solver.set_dc_map(dc_map);
    solver.setup_system(true);
    solver.assemble_system();
    solver.solve(0.5);

    ExpressionScalarPostprocessor<2> expression_postp_count_elem{"1", nu_map};
    ExpressionScalarPostprocessor<2> expression_postp_energy{"(Bx_q1^2+By_q1^2)/2.0 * nu_q1 * JxW_q1 + "
                                                             "(Bx_q2^2+By_q2^2)/2.0 * nu_q2 * JxW_q2 + "
                                                             "(Bx_q3^2+By_q3^2)/2.0 * nu_q3 * JxW_q3 + "
                                                             "(Bx_q4^2+By_q4^2)/2.0 * nu_q4 * JxW_q4", nu_map};

    double result_count_elem;
    expression_postp_count_elem.process(solver.get_triangulation(), solver.get_solution(), solver.get_fe(), result_count_elem);

    double result_energy;
    expression_postp_energy.process(solver.get_triangulation(), solver.get_solution(), solver.get_fe(), result_energy);

    std::cout << "result_count_elem: " << result_count_elem << std::endl;
    std::cout << "result_energy: " << result_energy << std::endl;

}

