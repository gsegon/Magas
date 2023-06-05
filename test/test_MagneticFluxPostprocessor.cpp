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
#include "MagneticFluxPostprocessor.h"

TEST(MagneticFluxPostprocessor, unit_square){

    std::string test_mesh = "/home/gordan/Programs/solver/test/test_data/test_unit_square/unit_square.msh";
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

    MagneticFluxPostprocessor<2> flux_postprocessor0{0, 0};
    MagneticFluxPostprocessor<2> flux_postprocessor1{0, 1};
    MagneticFluxPostprocessor<2> flux_postprocessor{0};
    std::vector<double> result;
    flux_postprocessor.process(solver.get_triangulation(), solver.get_solution(), solver.get_fe(), result);

    for (auto res: result)
        std::cout << res << ", " << std::endl;

}


