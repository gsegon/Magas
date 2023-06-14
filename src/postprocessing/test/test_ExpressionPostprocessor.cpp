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

    if((mat_id == 517), ( ((Bx_q1*By_q1*x_q1^2 +(By_q1^2-Bx_q1^2)*x_q1*y_q1-Bx_q1*By_q1*y_q1^2)/sqrt(x_q1^2+y_q1^2)*JxW_q1
                         +(Bx_q2*By_q2*x_q2^2 +(By_q2^2-Bx_q2^2)*x_q2*y_q2-Bx_q2*By_q2*y_q2^2)/sqrt(x_q2^2+y_q2^2)*JxW_q2
                         +(Bx_q3*By_q3*x_q3^2 +(By_q3^2-Bx_q3^2)*x_q3*y_q3-Bx_q3*By_q3*y_q3^2)/sqrt(x_q3^2+y_q3^2)*JxW_q3
                         +(Bx_q4*By_q4*x_q4^2 +(By_q4^2-Bx_q4^2)*x_q4*y_q4-Bx_q4*By_q4*y_q4^2)/sqrt(x_q4^2+y_q4^2)*JxW_q4) * nu/(46.4e-3-46.1e-3)), 0)

    ExpressionPostprocessor<2> expression_postp{"if(Bx_q1 >0.0, 1, 0)"};

    std::vector<double> result;
    expression_postp.process(solver.get_triangulation(), solver.get_solution(), solver.get_fe(), result);

    for (auto res: result)
        std::cout << res << ", " << std::endl;

}


