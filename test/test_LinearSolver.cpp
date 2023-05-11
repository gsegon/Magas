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
    solver.assemble_system();

}

TEST(LinearSolver, solve_system){

    LinearSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.setup_system();
    solver.assemble_system();
    solver.solve();

}

TEST(LinearSolver, output_results){

    LinearSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.setup_system();
    solver.assemble_system();
    solver.solve();
    solver.output_results("unit_square");

}
