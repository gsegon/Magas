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
#include <fstream>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/manifold_lib.h>
#include <fstream>
#include <deal.II/grid/grid_out.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/grid/grid_tools.h>


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

    constraints.clear();

    // Apply 0 DC boundary conditions
    for (auto& [mat_id, value] : dc_map){
        std::cout << "mat_id: " << mat_id << ": " << value << std::endl;
        VectorTools::interpolate_boundary_values(dof_handler, mat_id, Functions::ConstantFunction<2>(value), constraints);
    }

    // Apply periodic or anti-periodic boundary conditions
    for (auto [per_type, boundary_ids] : per_map){
        const IndexSet b_dofs_1 = DoFTools::extract_boundary_dofs(dof_handler, ComponentMask(), {boundary_ids[0]});
        const IndexSet b_dofs_2 = DoFTools::extract_boundary_dofs(dof_handler, ComponentMask(), {boundary_ids[1]});

        double weight = 0;
        if (per_type == "periodic")
            weight = 1;
        else if (per_type == "anti-periodic")
            weight = -1;

        AssertThrow ((per_type == "periodic") or (per_type == "anti-periodic"), ExceptionBase());
        AssertThrow (b_dofs_1.n_elements() == b_dofs_2.n_elements(), ExcInternalError());

        types::global_dof_index first;
        types::global_dof_index second;

        MappingQ1<dim> mapping;
        std::vector<Point<dim>> nodes(dof_handler.n_dofs());
        DoFTools::map_dofs_to_support_points(mapping, dof_handler, nodes);

        std::vector<std::pair<unsigned int, double>> dofs_1;
        for(auto dof : b_dofs_1){
            dofs_1.push_back({dof, nodes[dof].norm_square()});
        }

        std::vector<std::pair<unsigned int, double>> dofs_2;
        for(auto dof : b_dofs_2){
            dofs_2.push_back({dof, nodes[dof].norm_square()});
        }

        std::sort(dofs_1.begin(), dofs_1.end(), [](std::pair<unsigned int, double> a, std::pair<unsigned int, double> b) {return std::get<1>(a) < std::get<1>(b);});
        std::sort(dofs_2.begin(), dofs_2.end(), [](std::pair<unsigned int, double> a, std::pair<unsigned int, double> b) {return std::get<1>(a) < std::get<1>(b);});

        for (int i =0; i < dofs_1.size(); i++){

            first = std::get<0>(dofs_1[i]);
            second = std::get<0>(dofs_2[i]);

//            std::cout << "Setting dof " << first << " to dof " << second  << std::endl;
            if (abs(std::get<1>(dofs_1[i]) - std::get<1>(dofs_2[i])) > 1e-9)
                std::cerr << "Warning: Tol > 1e-8 " << std::endl;

            constraints.add_line(first);
            constraints.add_entry(first, second, weight);
        }

    }


//    constraints.print(std::cout);
//    std::ofstream dot_out("at_print.dot");
//    constraints.write_dot(dot_out);
    constraints.close();

    DynamicSparsityPattern dsp(dof_handler.n_dofs());

    // TODO: Investigate condensing DynamicSparsityPattern
//    constraints.condense(dsp);
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints);

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

        // Ass. local system into global
        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(cell_matrix,
                                               cell_rhs,
                                               local_dof_indices,
                                               system_matrix,
                                               system_rhs);
    }

    std::cout << "Number of dofs: " << dof_handler.n_dofs() << std::endl;

}

template<int dim>
void LinearSolver<dim>::solve(){

    // TODO: Investigate what happens here with constraints condense system_matrix and system_rhs.
    //  Whatever it is, it breaks periodicity
//    constraints.condense(system_matrix);
//    constraints.condense(system_rhs);

    SolverControl solver_control(10000, 1e-12);
    SolverCG<Vector<double>> solver(solver_control);

    PreconditionSSOR<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix, 1.6);

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
void LinearSolver<dim>::set_per_map(std::unordered_map<std::string, std::vector<unsigned int>> map) {
    this->per_map = map;
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

