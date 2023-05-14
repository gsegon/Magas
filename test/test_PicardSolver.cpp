//
// Created by gordan on 3/12/23.
//
#include <gtest/gtest.h>
#include <tuple>
#include <unordered_map>
#include <string>
#include <any>
#include <fstream>

#include "PicardSolver.h"
#include "PostMagneticFluxDensity.h"

TEST(PicardSolver, instantiation){

    PicardSolver<2> solver;

}

TEST(PicardSolver, read_mesh){

    std::string test_mesh = "/home/gordan/Programs/solver/test/test_data/test_EI_core/EI_core.msh";

    PicardSolver<2> solver;
    solver.read_mesh(test_mesh);


}

TEST(PicardSolver, setup_cell_nu_history){

    std::string test_mesh = "/home/gordan/Programs/solver/test/test_data/test_EI_core/EI_core.msh";

    PicardSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.setup_cell_nu_history();

}

TEST(PicardSolver, setup_system){

    std::string test_mesh = "/home/gordan/Programs/solver/test/test_data/test_EI_core/EI_core.msh";

    PicardSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.setup_cell_nu_history();
    solver.setup_system();

}

TEST(PicardSolver, assemble_system){

    std::string test_mesh = "/home/gordan/Programs/solver/test/test_data/test_unit_square/unit_square.msh";
    std::unordered_map<int, std::any> nu_map{{6, "Nonlinear"}};
    std::unordered_map<int, double> f_map{{6, 1.0}};
    std::unordered_map<int, double> dc_map{{5, 0.0}};

    PicardSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.setup_cell_nu_history();
    solver.setup_system();
    solver.set_nu_map(nu_map);
    solver.set_f_map(f_map);
    solver.set_dc_map(dc_map);
    solver.assemble_system();

}

TEST(PicardSolver, solve_system){

    std::string test_mesh = "/home/gordan/Programs/solver/test/test_data/test_unit_square/unit_square.msh";
    std::unordered_map<int, std::any> nu_map{{6, "Nonlinear"}};
    std::unordered_map<int, double> f_map{{6, 1.0}};
    std::unordered_map<int, double> dc_map{{5, 0.0}};

    PicardSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.setup_cell_nu_history();
    solver.setup_system();
    solver.set_nu_map(nu_map);
    solver.set_f_map(f_map);
    solver.set_dc_map(dc_map);
    solver.initialize_cell_nu_history(1);
    solver.assemble_system();
    solver.solve();

}

TEST(PicardSolver, output_results){

    std::string test_mesh = "/home/gordan/Programs/solver/test/test_data/test_unit_square/unit_square.msh";
    std::unordered_map<int, std::any> nu_map{{6, "Nonlinear"}};
    std::unordered_map<int, double> f_map{{6, 1.0}};
    std::unordered_map<int, double> dc_map{{5, 0.0}};

    PicardSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.setup_cell_nu_history();
    solver.setup_system();
    solver.set_nu_map(nu_map);
    solver.set_f_map(f_map);
    solver.set_dc_map(dc_map);
    solver.initialize_cell_nu_history(1);
    solver.assemble_system();
    solver.solve();

}

TEST(PicardSolver, solve_nonlinear){

    std::string test_mesh = "/home/gordan/Programs/solver/test/test_data/test_unit_square/unit_square.msh";
    std::unordered_map<int, std::any> nu_map{{6, "Nonlinear"}};
    std::unordered_map<int, double> f_map{{6, 1.0}};
    std::unordered_map<int, double> dc_map{{5, 0.0}};

    PicardSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.setup_cell_nu_history();
    solver.setup_system();
    solver.set_nu_map(nu_map);
    solver.set_f_map(f_map);
    solver.set_dc_map(dc_map);
    solver.initialize_cell_nu_history(1);
    solver.solve_nonlinear(3);

}


TEST(PicardSolver, EI_core){

    constexpr double mu_0 = 1.2566370621219e-6;
    constexpr double nu_0 = 1/mu_0;
    constexpr double nu_core = 1/(2500*mu_0);
    double J1 = 30*66/8.0645e-05;
    double J2 = -30*66/8.0645e-05;


    std::string test_mesh = "/home/gordan/Programs/solver/test/test_data/test_EI_core/EI_core.msh";

    std::unordered_map<int, std::any> nu_map{{200, "Nonlinear"},       // Core1
                                             {201, "Nonlinear"},       // Core2
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

    PicardSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.setup_cell_nu_history();
    solver.setup_system();
    solver.set_nu_map(nu_map);
    solver.set_f_map(f_map);
    solver.set_dc_map(dc_map);
    solver.initialize_cell_nu_history(nu_core);
    solver.solve_nonlinear(10);

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

    Vector<double> nu_val1(solver.get_triangulation().n_active_cells());
    Vector<double> nu_val2(solver.get_triangulation().n_active_cells());
    Vector<double> nu_val3(solver.get_triangulation().n_active_cells());
    Vector<double> nu_val4(solver.get_triangulation().n_active_cells());
    NuHistory<2> *local_quadrature_points_history;
    for (const auto& cell : dof_handler.active_cell_iterators()){

        if (nu_map.at(cell->material_id()).type() != typeid(double)){
            local_quadrature_points_history = reinterpret_cast<NuHistory<2> *>(cell->user_pointer());
            nu_val1(cell->active_cell_index()) = (double) local_quadrature_points_history->nu[0];
            nu_val2(cell->active_cell_index()) = (double) local_quadrature_points_history->nu[1];
            nu_val3(cell->active_cell_index()) = (double) local_quadrature_points_history->nu[2];
            nu_val4(cell->active_cell_index()) = (double) local_quadrature_points_history->nu[3];
        }
        else{
            nu_val1(cell->active_cell_index()) = 0*std::any_cast<double>(nu_map.at(cell->material_id()));
            nu_val2(cell->active_cell_index()) = 0*std::any_cast<double>(nu_map.at(cell->material_id()));
            nu_val3(cell->active_cell_index()) = 0*std::any_cast<double>(nu_map.at(cell->material_id()));
            nu_val4(cell->active_cell_index()) = 0*std::any_cast<double>(nu_map.at(cell->material_id()));
        }


    }
    data_out.add_data_vector(nu_val1, "nu_q1c", DataOut<2>::type_cell_data);
    data_out.add_data_vector(nu_val2, "nu_q2");
    data_out.add_data_vector(nu_val3, "nu_q3");
    data_out.add_data_vector(nu_val4, "nu_q4");

    data_out.build_patches();

    std::string filename = "test_result_EI_core_picard";
    std::ofstream output(filename + ".vtu");
    data_out.write_vtu(output);
}
