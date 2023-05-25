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
#include <deal.II/grid/tria_accessor.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/grid/grid_refinement.h>
#include <deal.II/numerics/error_estimator.h>

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

    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());

    constraints.clear();
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);

    for (auto& [mat_id, value] : dc_map){
        VectorTools::interpolate_boundary_values(dof_handler, mat_id, Functions::ConstantFunction<2>(value), constraints);
    }
    constraints.close();

    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler,
                                    dsp,
                                    constraints,
                                    false);

    sparsity_pattern.copy_from(dsp);
    system_matrix.reinit(sparsity_pattern);


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

        double nu = 0;
        double f = 0;
        Tensor<1, dim> Hc({0, 0});

        cell_matrix = 0.0;
        cell_rhs = 0.0;

        fe_values.reinit(cell);

        nu = nu_map.at(cell->material_id());
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
            for (const unsigned int i : fe_values.dof_indices()){
                for (const unsigned int j : fe_values.dof_indices()){
                    cell_matrix(i, j) += nu*                      // nu at cell
                            fe_values.shape_grad(i, q)*        // grad phi_i(x_q)
                            fe_values.shape_grad(j, q)*     // grad phi_j(x_q)
                            fe_values.JxW(q);                 // dx
                }

                // Current contribution
                cell_rhs(i) += f*                               // f at cell
                        fe_values.shape_value(i, q)*            // phi_i(x_q)
                        fe_values.JxW(q);

                // Magnet contribution
                cell_rhs(i) += Hc*fe_values.shape_grad(i, q)*fe_values.JxW(q);

            }
        }

        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
    }

    std::cout   << "\tNumber of degrees of freedom:\t"
                << dof_handler.n_dofs() << std::endl;
}

template <int dim>
void LinearSolver<dim>::refine_grid(){
    Vector<float> estimated_error_per_cell(triangulation.n_active_cells());

    KellyErrorEstimator<dim>::estimate(dof_handler,
                                       QGauss<dim -1>(fe.degree),
                                       {},
                                       solution,
                                       estimated_error_per_cell);


    GridRefinement::refine_and_coarsen_fixed_number(triangulation, estimated_error_per_cell, 0.05, 0.03);

    triangulation.execute_coarsening_and_refinement();
}

template<int dim>
void LinearSolver<dim>::solve(){

    SolverControl solver_control(10000, 1e-12);
    SolverCG<Vector<double>> solver(solver_control);

    PreconditionSSOR<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);

    solver.solve(system_matrix, solution, system_rhs, preconditioner);

    std::cout << "\t" << solver_control.last_step() << " CG iterations needed to obtain convergence." << std::endl;

    constraints.distribute(solution);

}

template<int dim>
void LinearSolver<dim>::set_nu_map(std::unordered_map<int, double> map) {
    this->nu_map = map;
}

template<int dim>
void LinearSolver<dim>::set_f_map(std::unordered_map<int, std::variant<double, std::pair<double, double>>> map) {
    this->f_map = map;
}

template<int dim>
void LinearSolver<dim>::set_dc_map(std::unordered_map<int, double> map) {
    this->dc_map = map;
}

template<int dim>
Vector<double>& LinearSolver<dim>::get_solution(){
    return this->solution;
}

template<int dim>
Vector<double>& LinearSolver<dim>::get_rhs(){
    return this->system_rhs;
}

template<int dim>
FE_Q<dim>& LinearSolver<dim>::get_fe(){
    return this->fe;
}

