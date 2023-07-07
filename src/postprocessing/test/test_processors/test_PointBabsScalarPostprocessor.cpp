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
#include "processors/PointBabsScalarPostprocessor.h"
#include "BHCurve.h"
#include "LinearBHCurve.h"

TEST(PointBabsScalarPostprocessor, torque_benchmark_kelvin_1){

    std::filesystem::path home = std::getenv("HOME");
    std::filesystem::path test_mesh = "../../../examples/torque_benchmark_kelvin_1/torque_benchmark_kelvin_1.msh";


    std::unordered_map<int, BHCurve*> nu_map{{3, new LinearBHCurve{795774.715025}},
                                           {4, new LinearBHCurve{795774.715025}},
                                           {5, new LinearBHCurve{795774.715025}},
                                           {6, new LinearBHCurve{795774.715025}},
                                           {8, new LinearBHCurve{795774.715025}}};

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

    PointBabsScalarPostprocessor<2> point_b_abs{{0, 0.003}};
    PointBabsScalarPostprocessor<2> point_b_abs2{{0, 15}};

    double result;
    point_b_abs.process(solver.get_triangulation(), solver.get_solution(), solver.get_fe(), result);
    std::cout << "result: " << result << std::endl;
    ASSERT_NEAR(result, 1.04108, 1e-3);


}




