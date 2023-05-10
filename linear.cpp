#include <iostream>
#include <deal.II/grid/tria.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/manifold_lib.h>
#include <fstream>

#include "linear.h"
using namespace dealii;

double RightHandSide::value(int physical_id) {
    return 0;
}

template<int dim>
LinearSolver<dim>::LinearSolver(){};

template<int dim>
void LinearSolver<dim>::read_mesh(std::string mesh_filepath) {
    GridIn<dim> grid_in;
    grid_in.attach_triangulation(triangulation);
    std::ifstream input_file(mesh_filepath);
    grid_in.read_msh(input_file);
    print_mesh_info(triangulation);

}

int main() {
    std::cout << "Hello, World!" << std::endl;
    return 0;
}
