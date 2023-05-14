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
    double J1 = 10*66/8.0645e-05;
    double J2 = -10*66/8.0645e-05;


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
    data_out.build_patches();

    std::string filename = "test_result_EI_core_picard";
    std::ofstream output(filename + ".vtu");
    data_out.write_vtu(output);
}
