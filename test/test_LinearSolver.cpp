//
// Created by gordan on 3/12/23.
//
#include <gtest/gtest.h>
#include <tuple>
#include <unordered_map>
#include <string>

#include "LinearSolver.h"

template<int dim>
void print_mesh_info(const Triangulation<dim> &triangulation)
{

    std::cout   << "Mesh info: " << std::endl
                << " dimension: " << dim << std::endl
                << " no. of cells: " << triangulation.n_active_cells() << std::endl;
    {
        std::map<types::boundary_id, unsigned int> boundary_count;
        for (const auto &face : triangulation.active_face_iterators())
            if (face->at_boundary())
                boundary_count[face->boundary_id()]++;

        std::cout << " boundary indicators: ";
        for (const std::pair<const types::boundary_id, unsigned int> &pair : boundary_count)
        {
            std::cout << pair.first << "(" << pair.second << " times)";
        }
        std::cout << std::endl;
    }
}

std::string test_mesh = "/home/gordan/Programs/solver/test/test_data/test_unit_square/unit_square.msh";
std::unordered_map<int, double> nu_map{{6, 1}};
std::unordered_map<int, double> f_map{{6, 1}};
std::unordered_map<int, double> dc_map{{5, 0}};


TEST(LinearSolver, instantiation){

    LinearSolver<2> solver;

}

TEST(LinearSolver, read_mesh){

    LinearSolver<2> solver;
    solver.read_mesh(test_mesh);
    print_mesh_info<2>(solver.get_triangulation());

}


TEST(LinearSolver, setup_system){

    LinearSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.setup_system();

}

TEST(LinearSolver, assemble_system){

    LinearSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.setup_system();
    solver.set_nu_map(nu_map);
    solver.set_f_map(f_map);
    solver.set_dc_map(dc_map);
    solver.assemble_system();

}

TEST(LinearSolver, solve_system){

    LinearSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.setup_system();
    solver.set_nu_map(nu_map);
    solver.set_f_map(f_map);
    solver.set_dc_map(dc_map);
    solver.assemble_system();
    solver.solve();

}

TEST(LinearSolver, output_results){

    LinearSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.setup_system();
    solver.set_nu_map(nu_map);
    solver.set_f_map(f_map);
    solver.set_dc_map(dc_map);
    solver.assemble_system();
    solver.solve();
    solver.output_results("unit_square");

}

TEST(LinearSolver, set_nu_map){

    LinearSolver<2> solver;
    solver.set_nu_map(nu_map);
}

TEST(LinearSolver, set_f_map){

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


    std::string test_mesh = "/home/gordan/Programs/solver/test/test_data/test_2_conductors/2_conductors_x_dense.msh";

    std::unordered_map<int, double> nu_map{{2, nu_0},
                                           {3, nu_0},
                                           {4, nu_0}};

    std::unordered_map<int, double> f_map{{4, 0},
                                          {2, Jdensity},
                                          {3, -Jdensity}};

    std::unordered_map<int, double> dc_map{{1, 0}};

    LinearSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.setup_system();
    solver.set_nu_map(nu_map);
    solver.set_f_map(f_map);
    solver.set_dc_map(dc_map);
    solver.assemble_system();
    solver.solve();
    solver.output_results("2_conductors");
    print_mesh_info(solver.get_triangulation());


}