// Magas - Magnetostatic Analysis Suite
// Copyright (C) 2023  Gordan Segon
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
// USA

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
#include <deal.II/grid/tria_accessor.h>
#include <fstream>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/function_cspline.h>

#include "include/PicardSolver.h"
using namespace dealii;

double h_fun_1(double b){
    double alpha = 2077.205761389225;
    double beta = 5.289952851132246;

    return alpha*b+std::exp(beta*b)-1;
}

double h_fun(double b){

    double nu_0 = 795774.715025;
    double b_sat = 2.2530727020352703;
    double theta = -1638220.518181392;

    if (b < b_sat)
        return h_fun_1(b);
    else
        return nu_0*b + theta;

}


template class PicardSolver<2>;

template<int dim>
PicardSolver<dim>::PicardSolver(): fe(1), dof_handler(triangulation), quadrature_formula(fe.degree+1)
{}

template<int dim>
void PicardSolver<dim>::read_mesh(const std::string& mesh_filepath) {
    GridIn<dim> grid_in;
    grid_in.attach_triangulation(triangulation);
    std::ifstream input_file(mesh_filepath);
    grid_in.read_msh(input_file);
}

template<int dim>
Triangulation<dim>& PicardSolver<dim>::get_triangulation(){
    return this->triangulation;
}

template<int dim>
void PicardSolver<dim>::setup_system() {
    dof_handler.distribute_dofs(fe);

    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);

}

template<int dim>
void PicardSolver<dim>::reinit_system(){
    system_matrix.reinit(sparsity_pattern);

    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
}

template<int dim>
void PicardSolver<dim>::assemble_system() {

    FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_gradients | update_quadrature_points |
                                                    update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double> cell_rhs(dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    double nu = 0;
    double f = 0;

    // Iterate over cells and assemble to local and move to global
    for (const auto &cell : dof_handler.active_cell_iterators()){

        cell_matrix = 0.0;
        cell_rhs = 0.0;

        fe_values.reinit(cell);

        // If nonlinear
        NuHistory<dim> *local_quadrature_points_history;
        if (nu_map.at(cell->material_id()).type() != typeid(double))
            local_quadrature_points_history = reinterpret_cast<NuHistory<dim> *>(cell->user_pointer());
        else
            nu = std::any_cast<double>(nu_map.at(cell->material_id()));

        f = f_map.at(cell->material_id());

        // Ass. into a local system
        for (const unsigned int q : fe_values.quadrature_point_indices()){

            if (nu_map.at(cell->material_id()).type() != typeid(double)){
                nu = local_quadrature_points_history->nu[q];
            }

            for (const unsigned int i : fe_values.dof_indices()){
                for (const unsigned int j : fe_values.dof_indices()){
                    cell_matrix(i, j) += nu*                      // nu at cell
                            fe_values.shape_grad(i, q)*        // grad phi_i(x_q)
                            fe_values.shape_grad(j, q)*     // grad phi_j(x_q)
                            fe_values.JxW(q);                 // dx
                }
                cell_rhs(i) += f*                               // f at cell
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
            system_rhs(local_dof_indices[i]) += cell_rhs(i);
        }
    }

    // Apply boundary conditions
    std::map<types::global_dof_index, double> boundary_values;
    for (auto& [mat_id, value] : dc_map){
        VectorTools::interpolate_boundary_values(dof_handler, mat_id, Functions::ConstantFunction<2>(value), boundary_values);
    }

    MatrixTools::apply_boundary_values(boundary_values,system_matrix,solution,system_rhs);

    std::cout << "Number of dofs: " << dof_handler.n_dofs() << std::endl;

}

template<int dim>
void PicardSolver<dim>::solve(){

    SolverControl solver_control(10000, 1e-12);
    SolverCG<Vector<double>> solver(solver_control);

    PreconditionSSOR<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix, 1.6);

    solver.solve(system_matrix, solution, system_rhs, preconditioner);

    std::cout << "\t" << solver_control.last_step() << " CG iterations needed to obtain convergence." << std::endl;

}

template<int dim>
void PicardSolver<dim>::solve_nonlinear(int max_iterations){
    while(max_iterations--){
        this->reinit_system();
        this->assemble_system();
        this->solve();
        this->update_cell_nu_history();
    }
}


template<int dim>
void PicardSolver<dim>::set_nu_map(std::unordered_map<int, std::any> map) {
    this->nu_map = map;
}

template<int dim>
void PicardSolver<dim>::set_f_map(std::unordered_map<int, double> map) {
    this->f_map = map;
}

template<int dim>
void PicardSolver<dim>::set_dc_map(std::unordered_map<int, double> map) {
    this->dc_map = map;
}

template<int dim>
Vector<double>& PicardSolver<dim>::get_solution(){
    return this->solution;
}

template<int dim>
Vector<double>& PicardSolver<dim>::get_rhs(){
    return this->system_rhs;
}

template<int dim>
FE_Q<dim>& PicardSolver<dim>::get_fe(){
    return this->fe;
}


template<int dim>
void PicardSolver<dim>::setup_cell_nu_history() {
    triangulation.clear_user_data();

    std::vector<NuHistory<dim>> tmp_storage;
    quadrature_nu_history.swap(tmp_storage);
    quadrature_nu_history.resize(triangulation.n_active_cells()*quadrature_formula.size());

    // Setup
    unsigned int history_index = 0;
    for (auto& cell : triangulation.active_cell_iterators()){
        cell->set_user_pointer(&quadrature_nu_history[history_index]);
        history_index += quadrature_formula.size();
    }
    Assert(history_index == quadrature_nu_history.size(), ExcInternalError());

}

template<int dim>
void PicardSolver<dim>::initialize_cell_nu_history(double initial_nu) {

    FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_gradients);
    for (auto& cell : triangulation.active_cell_iterators()){
        auto *local_quadrature_points_history = reinterpret_cast<NuHistory<dim> *>(cell->user_pointer());
        if (nu_map.at(cell->material_id()).type() != typeid(double))
            for (auto q: fe_values.quadrature_point_indices())
                local_quadrature_points_history->nu[q] = initial_nu;
    }
}


template<int dim>
void PicardSolver<dim>::update_cell_nu_history() {

    FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_gradients);
    std::vector<Tensor<1, dim>> solution_gradients(quadrature_formula.size());

    std::vector<double> b_at_qs(quadrature_formula.size(), 0);
    for (auto& cell : dof_handler.active_cell_iterators()){
        if (typeid(nu_map[cell->material_id()]) != typeid(double)){
            auto *local_quadrature_points_history = reinterpret_cast<NuHistory<dim> *>(cell->user_pointer());
            for (const unsigned int q: fe_values.quadrature_point_indices()){
                fe_values.reinit(cell);
                fe_values.get_function_gradients(solution, solution_gradients);
                b_at_qs[q] = std::sqrt(std::pow(solution_gradients[q][0], 2) + std::pow(solution_gradients[q][1], 2));

                double H_at_q = h_fun(b_at_qs[q]);
                double Nu_at_q = H_at_q/b_at_qs[q];

                local_quadrature_points_history->nu[q] += 0.1*(Nu_at_q-local_quadrature_points_history->nu[q]);
            }

        }
    }
}