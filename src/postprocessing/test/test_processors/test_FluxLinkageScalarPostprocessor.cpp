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
#include "processors/FluxLinkageScalarPostprocessor.h"
#include "NuCurve.h"
#include "LinearNuCurve.h"

TEST(FluxLinkageScalarPostprocessor, 6_conductors_3_circuits){

    std::string test_mesh = "../../../examples/6_conductors_3_circuits/6_conductors_3_circuits.msh";
    std::unordered_map<int, NuCurve*> nu_map{{2, new LinearNuCurve{795774.715025}},
                                             {3, new LinearNuCurve{795774.715025}},
                                             {4, new LinearNuCurve{795774.715025}},
                                             {5, new LinearNuCurve{795774.715025}},
                                             {6, new LinearNuCurve{795774.715025}},
                                             {7, new LinearNuCurve{795774.715025}},
                                             {8, new LinearNuCurve{795774.715025}}};

    std::unordered_map<int, std::variant<double, std::pair<double, double>>> f_map{{2, 0},
                                                                                   {3, 31.8309886184},
                                                                                   {4, -31.8309886184},
                                                                                   {5, -31.8309886184},
                                                                                   {6, 31.8309886184},
                                                                                   {7, -31.8309886184},
                                                                                   {8, 31.8309886184},
                                                                                   };
    std::unordered_map<int, double> dc_map{{1, 0}};

    LinearSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.set_nu_map(nu_map);
    solver.set_f_map(f_map);
    solver.set_dc_map(dc_map);
    solver.setup_system();
    solver.assemble_system();
    solver.solve();

    double result_4;
    FluxLinkageScalarPostprocessor<2> flux_linkage_pp_4{4, f_map};

    double result_8;
    FluxLinkageScalarPostprocessor<2> flux_linkage_pp_8{8, f_map};

    flux_linkage_pp_4.process(solver.get_triangulation(), solver.get_solution(), solver.get_fe(), result_4);
    flux_linkage_pp_8.process(solver.get_triangulation(), solver.get_solution(), solver.get_fe(), result_8);

    std::cout << "result_4: " << result_4 << std::endl;
    std::cout << "result_8: " << result_8 << std::endl;
    std::cout << "result_4 + result_8: " << result_4 + result_8 << std::endl;


}


