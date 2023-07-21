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
#include <deal.II/base/work_stream.h>

#include "NewtonSolver.h"
#include "PeriodicityMapperFactory.h"
#include "ConstFSource.h"
#include "Solver.h"

using namespace dealii;

template class NewtonSolver<2>;

template<int dim>
void NewtonSolver<dim>::setup_system(const bool initial_step) {

    if (initial_step) {
        Solver<dim>::dof_handler.distribute_dofs(Solver<dim>::fe);
        current_solution.reinit(Solver<dim>::dof_handler.n_dofs());

        Solver<dim>::constraints.clear();
        // Apply 0 DC boundary conditions
        for (auto& [mat_id, value] : Solver<dim>::dc_map){
            std::cout << "mat_id: " << mat_id << ": " << value << std::endl;
            VectorTools::interpolate_boundary_values(Solver<dim>::dof_handler, mat_id, Functions::ConstantFunction<2>(value), Solver<dim>::constraints);
        }

        // Apply periodic or anti-periodic boundary conditions
        for (auto [per_type, boundary_ids] : Solver<dim>::per_map){
            const IndexSet b_dofs_1 = DoFTools::extract_boundary_dofs(Solver<dim>::dof_handler, ComponentMask(), {boundary_ids[0]});
            const IndexSet b_dofs_2 = DoFTools::extract_boundary_dofs(Solver<dim>::dof_handler, ComponentMask(), {boundary_ids[1]});

            AssertThrow (b_dofs_1.n_elements() == b_dofs_2.n_elements(), ExcInternalError())

            std::vector<Point<dim>> nodes(Solver<dim>::dof_handler.n_dofs());
            DoFTools::map_dofs_to_support_points(MappingQ1<dim>(), Solver<dim>::dof_handler, nodes);

            std::map<unsigned int, std::vector<double>> dof_to_node;
            std::vector<unsigned int> dofs_1;
            std::vector<unsigned int> dofs_2;

            for (auto dof : b_dofs_1){
                dofs_1.push_back(dof);
                dof_to_node[dof] = {nodes[dof][0], nodes[dof][1]};
            }

            for (auto dof : b_dofs_2){
                dofs_2.push_back(dof);
                dof_to_node[dof] = {nodes[dof][0], nodes[dof][1]};
            }

            // Refactor to factory
            PeriodicityMapperFactory pmf{dofs_1, dofs_2, dof_to_node};
            auto pm = pmf.create(per_type);
            auto matched_pairs = pm->get_matched_pair_indices();

            for (auto matched_pair: matched_pairs){

                auto first = matched_pair.first;
                auto second = matched_pair.second;

                Solver<dim>::constraints.add_line(first);
                Solver<dim>::constraints.add_entry(first, second, pm->get_weigth());
            }
        }
        Solver<dim>::constraints.close();
    }

    DynamicSparsityPattern dsp(Solver<dim>::dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(Solver<dim>::dof_handler, dsp, Solver<dim>::constraints);
    Solver<dim>::sparsity_pattern.copy_from(dsp);

}

template<int dim>
void NewtonSolver<dim>::local_assemble_system(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                              typename Solver<dim>::AssemblyScratchData &scratch_data,
                                              typename Solver<dim>::AssemblyCopyData &copy_data) {

    double no = 0;
    FSource* f;
    ConstFSource f_zero{0};
    Tensor<1, dim> Hc({0, 0});

    const unsigned int dofs_per_cell = Solver<dim>::fe.n_dofs_per_cell();
    const unsigned int n_q_points = scratch_data.fe_values.get_quadrature().size();

    copy_data.cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
    copy_data.cell_rhs.reinit(dofs_per_cell);
    copy_data.local_dof_indices.resize(dofs_per_cell);

    std::vector<Tensor<1, dim>> old_solution_gradients(n_q_points);

    scratch_data.fe_values.reinit(cell);
    scratch_data.fe_values.get_function_gradients(current_solution, old_solution_gradients);

    auto f_variant = Solver<dim>::f_map.at(cell->material_id());
    if(std::holds_alternative<FSource*>(f_variant)){
        f = std::get<FSource*>(f_variant);
        Hc[0] = 0;
        Hc[1] = 0;
    }
    else if (std::holds_alternative<std::pair<double, double>>(f_variant)){
        f = &f_zero;
        Hc[0] = -std::get<std::pair<double, double>>(f_variant).second;
        Hc[1] = std::get<std::pair<double, double>>(f_variant).first;
    }
    else{
        throw std::runtime_error("Invalid source input for material_id=" + std::to_string(cell->material_id()) + ". Check problem definition file.");
    }

    for (unsigned int q = 0; q < n_q_points; q++){

        NuCurve* bh = Solver<dim>::nu_map.at(cell->material_id());
        double b_abs = std::sqrt(std::pow(old_solution_gradients[q][0],2) + std::pow(old_solution_gradients[q][1],2));
        no = bh->get_nu(b_abs) + bh->get_nu_prime(b_abs)*b_abs; // Newton::nu_fun(b_abs) + Newton::nu_fun_prime(b_abs)*b_abs;


        for (unsigned int i = 0; i < dofs_per_cell; i++){
            const auto &sd = scratch_data;
            for (unsigned int j = 0; j < dofs_per_cell; j++){
                copy_data.cell_matrix(i, j) += no*                      // nu at cell
                                               sd.fe_values.shape_grad(i, q)*        // grad phi_i(x_q)
                                               sd.fe_values.shape_grad(j, q)*     // grad phi_j(x_q)
                                               sd.fe_values.JxW(q);                 // dx
            }
            // Current + Magnet contribution - (...)
            // Current contribution
            auto vertex = cell->vertex(i);
            double f_val = f->get_value(vertex[0], vertex[1]);   // TODO: for all FSources fval += ...
            copy_data.cell_rhs(i) += ((f_val*sd.fe_values.shape_value(i, q)+Hc*sd.fe_values.shape_grad(i, q))                               // (f
                                      -                                   // -
                                      no*                                 // nu at cell*
                                      sd.fe_values.shape_grad(i, q)*       // phi_i(x_q)*
                                      old_solution_gradients[q]           //nabla u_n
                                     )*                                  // )*
                                     sd.fe_values.JxW(q);                   // dx

        }
    }

    cell->get_dof_indices(copy_data.local_dof_indices);
}

template<int dim>
void NewtonSolver<dim>::solve(const double alpha){

    SolverControl solver_control(10000, 1e-12);
    SolverCG<Vector<double>> solver(solver_control);

    PreconditionSSOR<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(Solver<dim>::system_matrix, 1.2);

    newton_update.reinit(Solver<dim>::dof_handler.n_dofs());
    solver.solve(Solver<dim>::system_matrix, newton_update, Solver<dim>::system_rhs, preconditioner);
    std::cout << "\t" << solver_control.last_step() << " CG iterations needed to obtain convergence." << std::endl;

    // Distribute constraints (periodic).
    Solver<dim>::constraints.distribute(newton_update);

    // Apply Newton step to current_solution
    current_solution.add(alpha, newton_update);

}

template <int dim>
double NewtonSolver<dim>::compute_residual(double alpha) const
{
    Vector<double> residual(Solver<dim>::dof_handler.n_dofs());

    Vector<double> evaluation_point(Solver<dim>::dof_handler.n_dofs());
    evaluation_point = current_solution;
    evaluation_point.add(alpha, newton_update);
    Solver<dim>::constraints.distribute(evaluation_point);

    FEValues<dim>     fe_values(Solver<dim>::fe,
                                Solver<dim>::quadrature_formula,
                                update_values | update_gradients | update_quadrature_points |
                                update_JxW_values);

    const unsigned int dofs_per_cell = Solver<dim>::fe.n_dofs_per_cell();
    const unsigned int n_q_points    = Solver<dim>::quadrature_formula.size();

    Vector<double>              cell_residual(dofs_per_cell);
    std::vector<Tensor<1, dim>> gradients(n_q_points);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    double no = 0;
    FSource* f;
    ConstFSource f_zero{0};
    Tensor<1, dim> Hc({0, 0});

    for (const auto &cell : Solver<dim>::dof_handler.active_cell_iterators()){
        cell_residual = 0;
        fe_values.reinit(cell);
        fe_values.get_function_gradients(evaluation_point, gradients);

        auto f_variant = Solver<dim>::f_map.at(cell->material_id());
        if(std::holds_alternative<FSource*>(f_variant)){
            f = std::get<FSource*>(f_variant);
            Hc[0] = 0;
            Hc[1] = 0;
        }
        else if (std::holds_alternative<std::pair<double, double>>(f_variant)){
            f = &f_zero;
            Hc[0] = -std::get<std::pair<double, double>>(f_variant).second;
            Hc[1] = std::get<std::pair<double, double>>(f_variant).first;
        }
        else{
            std::cout << "Something else?" << std::endl;
        }

        for (unsigned int q = 0; q < n_q_points; ++q){
            auto bh = Solver<dim>::nu_map.at(cell->material_id());
            double b_abs = std::sqrt(std::pow(gradients[q][0],2) + std::pow(gradients[q][1],2));
            no = bh->get_nu(b_abs) + bh->get_nu_prime(b_abs)*b_abs; // Newton::nu_fun(b_abs) + Newton::nu_fun_prime(b_abs)*b_abs;


            for (unsigned int i = 0; i < dofs_per_cell; ++i){
                // Current contribution
                auto vertex = cell->vertex(i);
                double f_val = f->get_value(vertex[0], vertex[1]);   // TODO: for all FSources fval += ...

                cell_residual(i) += ((f_val*fe_values.shape_value(i, q)+Hc*fe_values.shape_grad(i, q))
                                     - fe_values.shape_grad(i, q)         // \nabla \phi_i
                                       * no                               // * a_n
                                       * gradients[q])                    // * \nabla u_n
                                    * fe_values.JxW(q);                  // * dx
            }
        }

        cell->get_dof_indices(local_dof_indices);
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
            residual(local_dof_indices[i]) += cell_residual(i);
    }

    Solver<dim>::constraints.distribute(residual);

    for (types::global_dof_index i : DoFTools::extract_boundary_dofs(Solver<dim>::dof_handler))
        residual(i) = 0;

    return residual.l2_norm();
}

template<int dim>
void NewtonSolver<dim>::run(){

    std::cout << "::Running solver::" << std::endl;
    std::cout << "\tSetting up system...";
    setup_system(true);
    std::cout << "Done!" << std::endl;

    std::cout << "\tRunning Newton iteration...";
    double alpha=0.1;
    double res;
    double res_trial;
    std::vector<double> steps{1, 1/2.0, 1/4.0, 1/8.0, 1/16.0, 1/32.0};

    int i = 0;
    while(i++ < 100) {
        if (i > 1) {
            for (double alpha_trial: steps) {
                res_trial = compute_residual(alpha_trial);
                if (res_trial < res) {
                    alpha = alpha_trial;
                    break;
                }
            }
        }

        Solver<dim>::assemble_system();
        std::cout << "alpha = " << alpha << std::endl;
        solve(alpha);
        res = compute_residual(alpha);
        std::cout << "\tResidual(" << i << "): " << res << std::endl;

        if (res < 1e-6) {
            std::cout << "Converged!";
            break;
        }
        Solver<dim>::solution = current_solution;
        std::cout << "\tSolving: Done!" << std::endl;
    }
}