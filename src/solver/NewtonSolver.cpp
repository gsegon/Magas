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
#include <variant>

#include "NewtonSolver.h"
#include "AnalyticBHCurve.h"

using namespace dealii;

template class NewtonSolver<2>;

template<int dim>
NewtonSolver<dim>::NewtonSolver(): fe(1), dof_handler(triangulation), quadrature_formula(fe.degree + 1)
{
    bh = new AnalyticBHCurve{};
};

template<int dim>
void NewtonSolver<dim>::read_mesh(const std::string& mesh_filepath) {
    GridIn<dim> grid_in;
    grid_in.attach_triangulation(triangulation);
    std::ifstream input_file(mesh_filepath);
    grid_in.read_msh(input_file);
}

template<int dim>
Triangulation<dim>& NewtonSolver<dim>::get_triangulation(){
    return this->triangulation;
};

template<int dim>
void NewtonSolver<dim>::setup_system(const bool initial_step) {

    if (initial_step) {
        dof_handler.distribute_dofs(fe);
        current_solution.reinit(dof_handler.n_dofs());

        hanging_node_constraints.clear();
        DoFTools::make_hanging_node_constraints(dof_handler, hanging_node_constraints);

        hanging_node_constraints.close();
    }

    newton_update.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);

    system_matrix.reinit(sparsity_pattern);

}

template<int dim>
void NewtonSolver<dim>::assemble_system() {

    system_matrix = 0;
    system_rhs = 0;

    FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_gradients | update_quadrature_points |
                                                    update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned n_q_points = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double> cell_rhs(dofs_per_cell);

    std::vector<Tensor<1, dim>> old_solution_gradients(n_q_points);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    // Iterate over cells and assemble to local and move to global
    for (const auto &cell : dof_handler.active_cell_iterators()){

        double no = 0;
        double f = 0;
        Tensor<1, dim> Hc({0, 0});

        cell_matrix = 0.0;
        cell_rhs = 0.0;

        fe_values.reinit(cell);
        fe_values.get_function_gradients(current_solution, old_solution_gradients);

        auto f_variant = f_map.at(cell->material_id());
        if(std::holds_alternative<double>(f_variant)){
            f = std::get<double>(f_variant);
            Hc[0] = 0;
            Hc[1] = 0;
        }
        else if (std::holds_alternative<std::pair<double, double>>(f_variant)){
            f = 0;
            Hc[0] = -std::get<std::pair<double, double>>(f_variant).second;
            Hc[1] = std::get<std::pair<double, double>>(f_variant).first;
        }
        else{
            std::cout << "Something else?" << std::endl;
        }

        // Ass. into a local system
        for (const unsigned int q : fe_values.quadrature_point_indices()){

            if (nu_map.at(cell->material_id()).type() != typeid(double)){
                double b_abs = std::sqrt(std::pow(old_solution_gradients[q][0],2) + std::pow(old_solution_gradients[q][1],2));
                no = bh->get_nu(b_abs) + bh->get_nu_prime(b_abs)*b_abs; // Newton::nu_fun(b_abs) + Newton::nu_fun_prime(b_abs)*b_abs;
            }
            else
                no = std::any_cast<double>(nu_map.at(cell->material_id()));

            for (const unsigned int i : fe_values.dof_indices()){
                for (const unsigned int j : fe_values.dof_indices()){
                    cell_matrix(i, j) += no*                      // nu at cell
                                         fe_values.shape_grad(i, q)*        // grad phi_i(x_q)
                                         fe_values.shape_grad(j, q)*     // grad phi_j(x_q)
                                         fe_values.JxW(q);                 // dx
                }
                cell_rhs(i) += (f*fe_values.shape_value(i, q)                               // (f
                                -                                   // -
                                no*                                 // nu at cell*
                                fe_values.shape_grad(i, q)*       // phi_i(x_q)*
                                old_solution_gradients[q]           //nabla u_n
                               )*                                  // )*
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

    hanging_node_constraints.condense(system_matrix);
    hanging_node_constraints.condense(system_rhs);

    // Apply boundary conditions
    std::map<types::global_dof_index, double> boundary_values;
    for (auto& [mat_id, value] : dc_map){
        VectorTools::interpolate_boundary_values(dof_handler, mat_id, Functions::ConstantFunction<2>(value), boundary_values);
    }

    MatrixTools::apply_boundary_values(boundary_values,system_matrix,newton_update,system_rhs);

    std::cout << "Number of dofs: " << dof_handler.n_dofs() << std::endl;

}

template<int dim>
void NewtonSolver<dim>::solve(const double alpha){

    SolverControl solver_control(10000, 1e-12);
    SolverCG<Vector<double>> solver(solver_control);

    PreconditionSSOR<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);

    solver.solve(system_matrix, newton_update, system_rhs, preconditioner);
    std::cout << "\t" << solver_control.last_step() << " CG iterations needed to obtain convergence." << std::endl;

    hanging_node_constraints.distribute(newton_update);

    current_solution.add(alpha, newton_update);

}

template <int dim>
double NewtonSolver<dim>::compute_residual() const
{
    Vector<double> residual(dof_handler.n_dofs());

    Vector<double> evaluation_point(dof_handler.n_dofs());
    evaluation_point = current_solution;
    //evaluation_point.add(alpha, newton_update);

    FEValues<dim>     fe_values(fe,
                                quadrature_formula,
                                update_values | update_gradients | update_quadrature_points |
                                update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    Vector<double>              cell_residual(dofs_per_cell);
    std::vector<Tensor<1, dim>> gradients(n_q_points);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    double no = 0;
    double f = 0;
    Tensor<1, dim> Hc({0, 0});


    for (const auto &cell : dof_handler.active_cell_iterators()){
        cell_residual = 0;
        fe_values.reinit(cell);
        fe_values.get_function_gradients(evaluation_point, gradients);

        auto f_variant = f_map.at(cell->material_id());
        if(std::holds_alternative<double>(f_variant)){
            f = std::get<double>(f_variant);
            Hc[0] = 0;
            Hc[1] = 0;
        }
        else if (std::holds_alternative<std::pair<double, double>>(f_variant)){
            f = 0;
            Hc[0] = -std::get<std::pair<double, double>>(f_variant).second;
            Hc[1] = std::get<std::pair<double, double>>(f_variant).first;
        }
        else{
            std::cout << "Something else?" << std::endl;
        }


        for (unsigned int q = 0; q < n_q_points; ++q){
            if (nu_map.at(cell->material_id()).type() != typeid(double)){
                double b_abs = std::sqrt(std::pow(gradients[q][0],2) + std::pow(gradients[q][1],2));
                no = bh->get_nu(b_abs) + bh->get_nu_prime(b_abs)*b_abs;
            }
            else
                no = std::any_cast<double>(nu_map.at(cell->material_id()));

            for (unsigned int i = 0; i < dofs_per_cell; ++i)
                cell_residual(i) += (f*fe_values.shape_value(i, q)
                                     - fe_values.shape_grad(i, q)         // \nabla \phi_i
                                       * no                               // * a_n
                                       * gradients[q])                        // * \nabla u_n
                                    * fe_values.JxW(q);                // * dx
        }

        cell->get_dof_indices(local_dof_indices);
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
            residual(local_dof_indices[i]) += cell_residual(i);
    }

    hanging_node_constraints.condense(residual);

    for (types::global_dof_index i : DoFTools::extract_boundary_dofs(dof_handler))
        residual(i) = 0;

    return residual.l2_norm();
}


//template<int dim>
//void NewtonSolver<dim>::solve_nonlinear(int max_iterations){
//
//    while(max_iterations--){
//        assemble_system();
//        double last_residual_norm = system_rhs.l2_norm();
//        solve();
//
//        std::cout << "  Residual: " << compute_residual(0) << std::endl;
//    }
//}


template<int dim>
void NewtonSolver<dim>::set_nu_map(std::unordered_map<int, std::any> map) {
    this->nu_map = map;
}

template<int dim>
void NewtonSolver<dim>::set_f_map(std::unordered_map<int, std::variant<double, std::pair<double, double>>> map) {
    this->f_map = map;
}

template<int dim>
void NewtonSolver<dim>::set_dc_map(std::unordered_map<int, double> map) {
    this->dc_map = map;
}

template<int dim>
Vector<double>& NewtonSolver<dim>::get_solution(){
    return this->solution;
}

template<int dim>
Vector<double>& NewtonSolver<dim>::get_current_solution(){
    return this->current_solution;
}

template<int dim>
Vector<double>& NewtonSolver<dim>::get_rhs(){
    return this->system_rhs;
}

template<int dim>
FE_Q<dim>& NewtonSolver<dim>::get_fe(){
    return this->fe;
}
