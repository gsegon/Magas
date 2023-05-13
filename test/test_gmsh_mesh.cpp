//
// Created by gordan on 3/12/23.
//
#include <gtest/gtest.h>
#include <tuple>
#include <unordered_map>
#include <string>
#include <fstream>

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
#include <deal.II/grid/manifold_lib.h>

using namespace dealii;

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

        std::cout << " boundary indicator: " << std::endl;
        for (const std::pair<const types::boundary_id, unsigned int> &pair : boundary_count)
        {
            std::cout << "\t\t" << pair.first << " (" << pair.second << " times)" << std::endl;
        }
        std::cout << std::endl;
    }
}


TEST(GmshMesh, unit_square){

    std::string test_mesh = "/home/gordan/Programs/solver/test/test_data/test_gmsh_mesh/unit_square.msh";

    Triangulation<2> triangulation;
    GridIn<2> grid_in;
    grid_in.attach_triangulation(triangulation);
    std::ifstream input_file(test_mesh);
    grid_in.read_msh(input_file);
    print_mesh_info(triangulation);

}

TEST(GmshMesh, two_conductors){

    std::string test_mesh = "/home/gordan/Programs/solver/test/test_data/test_gmsh_mesh/2_conductors_x.msh";

    Triangulation<2> triangulation;
    GridIn<2> grid_in;
    grid_in.attach_triangulation(triangulation);
    std::ifstream input_file(test_mesh);
    grid_in.read_msh(input_file);
    print_mesh_info(triangulation);


}

TEST(GmshMesh, test_n){

    std::string test_mesh = "/home/gordan/Programs/solver/test/test_data/test_gmsh_mesh/test_n.msh";

    Triangulation<2> triangulation;
    GridIn<2> grid_in;
    grid_in.attach_triangulation(triangulation);
    std::ifstream input_file(test_mesh);
    grid_in.read_msh(input_file);
    print_mesh_info(triangulation);

}

