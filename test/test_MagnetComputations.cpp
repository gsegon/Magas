//
// Created by gordan on 3/12/23.
//
#include <gtest/gtest.h>
#include <tuple>
#include <unordered_map>
#include <string>
#include <any>
#include <fstream>

#include <deal.II/grid/tria.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>

using namespace dealii;

TEST(MagnetComputations, concepts){

    std::string test_mesh = "/home/gordan/Programs/solver/test/test_data/test_magnet/BlockMagnet.msh";

    Triangulation<2> triangulation;
    GridIn<2> grid_in;
    grid_in.attach_triangulation(triangulation);
    std::ifstream input_file(test_mesh);
    grid_in.read_msh(input_file);

    FE_Q<2> fe(1);
    DoFHandler<2> dof_handler(triangulation);
    QGauss<2> quadrature_formula(fe.degree+1);

    dof_handler.distribute_dofs(fe);
    FEValues<2> fe_values(fe, quadrature_formula, update_values | update_gradients | update_quadrature_points |
                                                    update_JxW_values);

    int count = 0;
    for (auto cell : dof_handler.active_cell_iterators()){
        std::cout << "Number of lines: " << cell->n_lines() << std::endl;
        std::cout << "Lines: " << std::endl;
        for (auto line_index : cell->line_indices()){
            std::cout << line_index << ": ";
            auto line = cell->line(line_index);
            for (auto vertex_index : line->vertex_indices()){
                    auto vertex = line->vertex(vertex_index);
                    std::cout << vertex << ", ";
                }

            std::cout << "measure: " << line->measure() << std::endl;
            }

        std::cout << std::endl;
        count++;
        if(count > 3)
            break;
        }



}

