#include <iostream>
#include <fstream>

#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/lac/sparse_direct.h>

#include "PeriodicityMapperFactory.h"
#include "LinearSolver.h"
#include "ConstFSource.h"

using namespace dealii;

template class LinearSolver<2>;

template<int dim>
void LinearSolver<dim>::setup_system() {
    this->dof_handler.distribute_dofs(this->fe);

    this->constraints.clear();

    // setup_rotation
    for (auto [key, value] : this->rot_map){
        this->setup_rotation(key.first, key.second, value);
    }

    // Apply 0 DC boundary conditions
    for (auto& [mat_id, value] : this->dc_map){
        VectorTools::interpolate_boundary_values(this->dof_handler, mat_id, Functions::ConstantFunction<2>(value), this->constraints);
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

            if (first != second){
                Solver<dim>::constraints.add_line(first);
                Solver<dim>::constraints.add_entry(first, second, pm->get_weigth());
            }
        }
    }

//    constraints.print(std::cout);
//    std::ofstream dot_out("at_print.dot");
//    constraints.write_dot(dot_out);
    Solver<dim>::constraints.close();

    DynamicSparsityPattern dsp(Solver<dim>::dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(Solver<dim>::dof_handler, dsp, Solver<dim>::constraints);

    //Extend dsp due to rotation mappings:
    this->extend_dsp(dsp);

    Solver<dim>::sparsity_pattern.copy_from(dsp);

    Solver<dim>::system_matrix.reinit(Solver<dim>::sparsity_pattern);
    Solver<dim>::solution.reinit(Solver<dim>::dof_handler.n_dofs());
    Solver<dim>::system_rhs.reinit(Solver<dim>::dof_handler.n_dofs());
}

template<int dim>
void LinearSolver<dim>::local_assemble_system(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                              typename Solver<dim>::AssemblyScratchData &scratch_data,
                                              typename Solver<dim>::AssemblyCopyData &copy_data) {

    double nu = 0;
    FSource* f;
    ConstFSource f_zero{0};
    Tensor<1, dim> Hc({0, 0});

    const unsigned int dofs_per_cell = this->fe.n_dofs_per_cell();
    const unsigned int n_q_points = scratch_data.fe_values.get_quadrature().size();

    copy_data.cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
    copy_data.cell_rhs.reinit(dofs_per_cell);
    copy_data.local_dof_indices.resize(dofs_per_cell);

    scratch_data.fe_values.reinit(cell);

    NuCurve* bh = this->nu_map.at(cell->material_id());
    nu = bh->get_nu(0);
    auto f_variant = this->f_map.at(cell->material_id());
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

    for (unsigned int q = 0; q < n_q_points; q++)
        for (unsigned int i = 0; i < dofs_per_cell; i++){
            const auto &sd = scratch_data;
            for (unsigned int j = 0; j < dofs_per_cell; j++){
                copy_data.cell_matrix(i, j) += nu*              // nu at cell
                            sd.fe_values.shape_grad(i, q)*      // grad phi_i(x_q)
                            sd.fe_values.shape_grad(j, q)*     // grad phi_j(x_q)
                            sd.fe_values.JxW(q);               // dx
            }
            // Current contribution
            auto vertex = cell->vertex(i);
            double f_val = f->get_value(vertex[0], vertex[1]);   // TODO: for all FSources fval += ...
            copy_data.cell_rhs(i) += f_val*                     // f at cell
                    sd.fe_values.shape_value(i, q)*             // phi_i(x_q)
                    sd.fe_values.JxW(q);                        // dx

            // Magnet contribution
                copy_data.cell_rhs(i) += Hc*sd.fe_values.shape_grad(i, q)*sd.fe_values.JxW(q);
        }

    cell->get_dof_indices(copy_data.local_dof_indices);
}

template<int dim>
void LinearSolver<dim>::solve(){

    SparseDirectUMFPACK A_direct;
    this->solution = this->system_rhs;
    A_direct.solve(this->system_matrix, this->solution);

    this->constraints.distribute(this->solution);

}

template<int dim>
void LinearSolver<dim>::run(){
    std::cout << "::Running solver::" << std::endl;
    std::cout << "\tSetting up system...";
    setup_system();
    std::cout << "Done!" << std::endl;
    std::cout << "\tAssembling system...";
    this->assemble_system();
    std::cout << "Done!" << std::endl;
    std::cout << "\tSolving system...";
    solve();
    std::cout << "\tSolving: Done!" << std::endl;
}

