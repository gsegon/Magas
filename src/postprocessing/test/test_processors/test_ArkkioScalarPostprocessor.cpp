//
// Created by gordan on 3/12/23.
//
#include <gtest/gtest.h>
#include <tuple>
#include <unordered_map>
#include <string>
#include <any>
#include <fstream>

#include "LinearSolver.h"
#include "processors/ArkkioScalarPostprocessor.h"

TEST(ArkkioScalarPostprocessor, torque_benchmark_kelvin_1){

    std::string test_mesh = "~/Programs/solver/examples/torque_benchmark_kelvin_1/torque_benchmark_kelvin_1.msh";
    std::unordered_map<int, double> nu_map{{3, 795774.715025},
                                           {4, 795774.715025},
                                           {5, 795774.715025},
                                           {6, 795774.715025},
                                           {8, 795774.715025}};

    std::unordered_map<int, std::variant<double, std::pair<double, double>>> f_map{{3, std::pair<double, double>{0.0, 1000000.0}},
                                                                                   {4, 0},
                                                                                   {5, 0},
                                                                                   {6, 0},
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

    ArkkioScalarPostprocessor<2> arkkio{5, nu_map};

    double result;
    arkkio.process(solver.get_triangulation(), solver.get_solution(), solver.get_fe(), result);

    std::cout << "result: " << result << std::endl;
    std::cout << "result*2e-2: " << result*2e-2 << std::endl;

    ASSERT_NEAR(result*2e-2, 1.0, 1e-3);

}


