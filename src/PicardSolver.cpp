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

#include "PicardSolver.h"
using namespace dealii;

std::vector<double> B = {0,        0.2003,     0.3204,
                         0.40045,  0.50055,    0.5606,
                         0.7908,   0.931,      1.1014,
                         1.2016,   1.302,      1.4028,
                         1.524,    1.626,      1.698,
                         1.73,     1.87,       1.99,
                         2.04,     2.07,       2.095,
                         2.2,      2.4};
std::vector<double> H = {0,        238.7,     318.3,
                         358.1,    437.7,      477.5,
                         636.6,    795.8,      1114.1,
                         1273.2,   1591.5,     2228.2,
                         3183.1,   4774.6,     6366.2,
                         7957.7,   15915.5,    31831,
                         47746.5,  63663,      79577.5,
                         159155,   318310};

dealii::Functions::CSpline<1> sp(B, H);


template class PicardSolver<2>;

template<int dim>
PicardSolver<dim>::PicardSolver(): fe(1), dof_handler(triangulation), quadrature_formula(fe.degree+1)
{};

template<int dim>
void PicardSolver<dim>::read_mesh(std::string mesh_filepath) {
    GridIn<dim> grid_in;
    grid_in.attach_triangulation(triangulation);
    std::ifstream input_file(mesh_filepath);
    grid_in.read_msh(input_file);
}

template<int dim>
Triangulation<dim>& PicardSolver<dim>::get_triangulation(){
    return this->triangulation;
};

template<int dim>
void PicardSolver<dim>::setup_system() {
    dof_handler.distribute_dofs(fe);

    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);

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

    SolverControl solver_control(1000, 1e-12);
    SolverCG<Vector<double>> solver(solver_control);

    PreconditionSSOR<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix, 1.6);

    solver.solve(system_matrix, solution, system_rhs, preconditioner);

    std::cout << "\t" << solver_control.last_step() << " CG iterations needed to obtain convergence." << std::endl;

}

template<int dim>
void PicardSolver<dim>::solve_nonlinear(int max_iterations){
    while(max_iterations--){
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

    double Bmax = 2.3;
    for (auto& cell : dof_handler.active_cell_iterators()){
        if (typeid(nu_map[cell->material_id()]) != typeid(double)){
            for (const unsigned int q: fe_values.quadrature_point_indices()){
                auto *local_quadrature_points_history = reinterpret_cast<NuHistory<dim> *>(cell->user_pointer());
                fe_values.reinit(cell);
                fe_values.get_function_gradients(solution, solution_gradients);
                double Blocal = std::sqrt(std::pow(solution_gradients[q][0], 2) + std::pow(solution_gradients[q][1], 2));

                if (Blocal > Bmax)
                    Blocal = Bmax;

                double nu_vcalc = sp.value(dealii::Point<1>{Blocal});

                local_quadrature_points_history->nu[q] += 0.7*(nu_vcalc-local_quadrature_points_history->nu[q]);
            }
        }
    }
}