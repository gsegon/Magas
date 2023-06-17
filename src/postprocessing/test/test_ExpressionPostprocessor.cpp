//
// Created by gordan on 3/12/23.
//
#include <gtest/gtest.h>
#include <tuple>
#include <unordered_map>
#include <string>
#include <any>
#include <fstream>

#include "../../solver/include/LinearSolver.h"
#include "ExpressionPostprocessor.h"

TEST(ExpressionPostprocessor, unit_square){

    std::string test_mesh = "/home/gordan/Programs/solver/examples/unit_square/unit_square.msh";
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

    ExpressionPostprocessor<2> expression_postp{"if(Bx_q1 >0.0, 1, 0)"};

    std::vector<double> result;
    expression_postp.process(solver.get_triangulation(), solver.get_solution(), solver.get_fe(), result);

    for (auto res: result)
        std::cout << res << ", " << std::endl;

}


