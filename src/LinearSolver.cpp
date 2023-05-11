#include <iostream>
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
#include <deal.II/grid/manifold_lib.h>
#include <fstream>

#include "LinearSolver.h"
using namespace dealii;


template class LinearSolver<2>;

template<int dim>
LinearSolver<dim>::LinearSolver(): fe(1), dof_handler(triangulation), quadrature_formula(fe.degree+1)
{};

template<int dim>
void LinearSolver<dim>::read_mesh(std::string mesh_filepath) {
    GridIn<dim> grid_in;
    grid_in.attach_triangulation(triangulation);
    std::ifstream input_file(mesh_filepath);
    grid_in.read_msh(input_file);
}

template<int dim>
Triangulation<dim>& LinearSolver<dim>::get_triangulation(){
    return this->triangulation;
};

template<int dim>
void LinearSolver<dim>::setup_system() {
    dof_handler.distribute_dofs(fe);

    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);

    system_matrix.reinit(sparsity_pattern);

    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());

}

template<int dim>
void LinearSolver<dim>::assemble_system() {

    FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_gradients | update_quadrature_points |
                                                    update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double> cell_rhs(dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    // Iterate over cells and assemble to local and move to global
    for (const auto &cell : dof_handler.active_cell_iterators()){

        cell_matrix = 0.0;
        cell_rhs = 0.0;

        fe_values.reinit(cell);

        // Ass. into a local system
        for (const unsigned int q : fe_values.quadrature_point_indices()){
            for (const unsigned int i : fe_values.dof_indices()){
                for (const unsigned int j : fe_values.dof_indices()){
                    cell_matrix(i, j) +=
                            fe_values.shape_grad(i, q)*        // grad phi_i(x_q)
                            fe_values.shape_grad(j, q)*     // grad phi_j(x_q)
                            fe_values.JxW(q);                 // dx
                }
                cell_rhs(i) +=
                        fe_values.shape_value(i, q)*        // phi_i(x_q)
                        fe_values.JxW(q);                   // dx
            }
        }

        // Ass. local system into global
        cell->get_dof_indices(local_dof_indices);
        for (const unsigned int i : fe_values.dof_indices()){
            for (const unsigned int j : fe_values.dof_indices()){
                system_matrix.add(local_dof_indices[i], local_dof_indices[j], cell_matrix(i, j));
            }
            system_rhs.add(local_dof_indices[i]);
        }
    }

    // Apply boundary conditions
    std::map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values(dof_handler, 5, Functions::ZeroFunction<2>(), boundary_values);
    MatrixTools::apply_boundary_values(boundary_values,system_matrix,solution,system_rhs);

}

template<int dim>
void LinearSolver<dim>::solve(){

    SolverControl solver_control(1000, 1e-12);
    SolverCG<Vector<double>> solver(solver_control);

    PreconditionSSOR<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);

    solver.solve(system_matrix, solution, system_rhs, preconditioner);

    std::cout << "\t" << solver_control.last_step() << " CG iterations needed to obtain convergence." << std::endl;

}