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
#include "PostMagneticFluxDensity.h"

TEST(PostMagneticFluxDensity, initialize){

    std::string test_mesh = "/home/gordan/Programs/solver/test/test_data/test_unit_square/unit_square.msh";
    std::unordered_map<int, double> nu_map{{6, 1}};
    std::unordered_map<int, double> f_map{{6, 1}};
    std::unordered_map<int, double> dc_map{{5, 0}};

    LinearSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.setup_system();
    solver.set_nu_map(nu_map);
    solver.set_f_map(f_map);
    solver.set_dc_map(dc_map);
    solver.assemble_system();
    solver.solve();

    PostMagneticFluxDensity<2> pmfdp;

    DoFHandler<2> dof_handler(solver.get_triangulation());
    dof_handler.distribute_dofs(solver.get_fe());

    DataOut<2> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solver.get_solution(), "u");
    data_out.add_data_vector(solver.get_solution(), pmfdp);
    data_out.build_patches();

    std::string filename = "pmfdp";
    std::ofstream output(filename + ".vtu");
    data_out.write_vtu(output);

}

