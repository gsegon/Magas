//
// Created by gordan on 3/12/23.
//
#include <gtest/gtest.h>
#include <tuple>
#include <unordered_map>
#include <string>
#include <fstream>

#include "LinearSolver.h"
#include "PostMagneticFluxDensity.h"
#include <deal.II/numerics/data_out.h>

using namespace dealii;

TEST(LinearSolver, cell_out){

    constexpr double mu_0 = 1.2566370621219e-6;
    constexpr double nu_0 = 1/mu_0;
    constexpr double nu_core = 1/(2500*mu_0);
    double J1 = 10*66/8.0645e-05;
    double J2 = -10*66/8.0645e-05;


    std::string test_mesh = "/home/gordan/Programs/solver/test/test_data/test_EI_core/EI_core.msh";
    std::unordered_map<int, double> nu_map{{200, nu_core},       // Core1
                                           {201, nu_core},       // Core2
                                           {202, nu_0},       // Copper
                                           {203, nu_0},       // Copper
                                           {204, nu_0},       // Air
                                           {205, nu_0},       // Air
    };

    std::unordered_map<int, double> f_map{ {200, 0},        // Core1
                                           {201, 0},        // Core2
                                           {202, J1},       // Copper
                                           {203, J2},       // Copper
                                           {204, 0},        // Air
                                           {205, 0},        // Air
    };

    std::unordered_map<int, double> dc_map{{44, 0}} ;      // Air

    // Solve
    LinearSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.setup_system();
    solver.set_nu_map(nu_map);
    solver.set_f_map(f_map);
    solver.set_dc_map(dc_map);
    solver.assemble_system();
    solver.solve();

    // Post process
    PostMagneticFluxDensity<2> pmfdp;
    DoFHandler<2> dof_handler(solver.get_triangulation());
    dof_handler.distribute_dofs(solver.get_fe());
    DataOut<2> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solver.get_solution(), "u");
    data_out.add_data_vector(solver.get_solution(), pmfdp);

    Vector<float> mat_id_mask(solver.get_triangulation().n_active_cells());
    for (const auto& cell : dof_handler.active_cell_iterators()){
        mat_id_mask(cell->active_cell_index()) = (float) cell->material_id();
    }
    data_out.add_data_vector(mat_id_mask, "mat_id");

    Vector<float> at_boundary(solver.get_triangulation().n_active_cells());
    for (const auto& cell : dof_handler.active_cell_iterators()){
        at_boundary(cell->active_cell_index()) = (float) cell->at_boundary();
    }
    data_out.add_data_vector(at_boundary, "at_boundary");

    Vector<float> manifold_id(solver.get_triangulation().n_active_cells());
    for (const auto& cell : dof_handler.active_cell_iterators()){
        manifold_id(cell->active_cell_index()) = (float) cell->manifold_id();
    }
    data_out.add_data_vector(manifold_id, "manifold_id");

    data_out.build_patches();

    // Write to file
    std::string filename = "test_result_cell_out";
    std::ofstream output(filename + ".vtu");
    data_out.write_vtu(output);

//    std::ofstream output_(filename + ".vtu");
//    data_out.write_vtu(output);
}



TEST(LinearSolver, instantiation){

    LinearSolver<2> solver;

}

TEST(LinearSolver, read_mesh){
    std::string test_mesh = "/home/gordan/Programs/solver/test/test_data/test_unit_square/unit_square.msh";

    LinearSolver<2> solver;
    solver.read_mesh(test_mesh);


}

TEST(LinearSolver, setup_system){

    std::string test_mesh = "/home/gordan/Programs/solver/test/test_data/test_unit_square/unit_square.msh";

    LinearSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.setup_system();

}

TEST(LinearSolver, assemble_system){

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

}

TEST(LinearSolver, solve_system){

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

}

TEST(LinearSolver, set_nu_map){

    std::string test_mesh = "/home/gordan/Programs/solver/test/test_data/test_unit_square/unit_square.msh";
    std::unordered_map<int, double> nu_map{{6, 1}};
    std::unordered_map<int, double> f_map{{6, 1}};
    std::unordered_map<int, double> dc_map{{5, 0}};

    LinearSolver<2> solver;
    solver.set_nu_map(nu_map);
}

TEST(LinearSolver, set_f_map){

    std::string test_mesh = "/home/gordan/Programs/solver/test/test_data/test_unit_square/unit_square.msh";
    std::unordered_map<int, double> nu_map{{6, 1}};
    std::unordered_map<int, double> f_map{{6, 1}};
    std::unordered_map<int, double> dc_map{{5, 0}};

    LinearSolver<2> solver;
    solver.set_f_map(f_map);
}

TEST(LinearSolver, setup_system_3){

    std::string test_mesh = "/home/gordan/Programs/solver/test/test_data/test_2_conductors/2_conductors_x.msh";

    LinearSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.setup_system();

}

TEST(LinearSolver, 2_conductors){

    constexpr double mu_0 = 1.2566370621219e-6;
    constexpr double nu_0 = 1/mu_0;
    double i_current = 1e3;
    double Jdensity = i_current / (std::pow(0.1,2) * M_PI);


    std::string test_mesh = "/home/gordan/Programs/solver/test/test_data/test_2_conductors/2_conductors_x.msh";

    std::unordered_map<int, double> nu_map{{1, nu_0},       // Conductor 1
                                           {2, nu_0},       // Conductor 2
                                           {3, nu_0}};      // Air          

    std::unordered_map<int, double> f_map{ {1, Jdensity},   // Conductor 1
                                           {2, -Jdensity},  // Conductor 2
                                           {3, 0},          // Air
                                        };

    std::unordered_map<int, double> dc_map{{100, 0}};       // Outerbounds

    LinearSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.setup_system();
    solver.set_nu_map(nu_map);
    solver.set_f_map(f_map);
    solver.set_dc_map(dc_map);
    solver.assemble_system();
    solver.solve();

}

TEST(LinearSolver, EI_core){

    constexpr double mu_0 = 1.2566370621219e-6;
    constexpr double nu_0 = 1/mu_0;
    constexpr double nu_core = 1/(2500*mu_0);
    double J1 = 10*66/8.0645e-05;
    double J2 = -10*66/8.0645e-05;


    std::string test_mesh = "/home/gordan/Programs/solver/test/test_data/test_EI_core/EI_core.msh";
    std::unordered_map<int, double> nu_map{{200, nu_core},       // Core1
                                             {201, nu_core},       // Core2
                                             {202, nu_0},       // Copper
                                             {203, nu_0},       // Copper
                                             {204, nu_0},       // Air
                                             {205, nu_0},       // Air
    };

    std::unordered_map<int, double> f_map{ {200, 0},        // Core1
                                           {201, 0},        // Core2
                                           {202, J1},       // Copper
                                           {203, J2},       // Copper
                                           {204, 0},        // Air
                                           {205, 0},        // Air
    };

    std::unordered_map<int, double> dc_map{{44, 0}} ;      // Air

    // Solve
    LinearSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.setup_system();
    solver.set_nu_map(nu_map);
    solver.set_f_map(f_map);
    solver.set_dc_map(dc_map);
    solver.assemble_system();
    solver.solve();

    // Post process
    PostMagneticFluxDensity<2> pmfdp;
    DoFHandler<2> dof_handler(solver.get_triangulation());
    dof_handler.distribute_dofs(solver.get_fe());
    DataOut<2> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solver.get_solution(), "u");
    data_out.add_data_vector(solver.get_solution(), pmfdp);
    data_out.build_patches();

    // Write to file
    std::string filename = "test_result_EI_core_linear";
    std::ofstream output(filename + ".vtu");
    data_out.write_vtu(output);
}

