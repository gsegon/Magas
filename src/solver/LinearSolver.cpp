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

#include <deal.II/base/work_stream.h>
#include <deal.II/base/multithread_info.h>

#include "LinearSolver.h"
//#include "CirclePeriodicityMapper.h"
#include "LinePeriodicityMapper.h"
#include "CirclePeriodicityMapper.h"
#include "IPeriodicityMapper.h"
//#include "OverlapPointsTransformation.h"

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

        //TODO: Move periodic mapping outside, in a separate class/utility.
        double weight = 0;
        if (per_type == "periodic-line" or per_type == "periodic-circle")
            weight = 1;
        else if (per_type == "anti-periodic-line" or per_type == "anti-periodic-circle")
            weight = -1;

        AssertThrow (   (per_type == "periodic-line") or
                        (per_type == "anti-periodic-line") or
                        (per_type == "periodic-circle") or
                        (per_type == "anti-periodic-circle"), ExceptionBase())
        AssertThrow (b_dofs_1.n_elements() == b_dofs_2.n_elements(), ExcInternalError())

        std::vector<Point<dim>> nodes(dof_handler.n_dofs());
        DoFTools::map_dofs_to_support_points(MappingQ1<dim>(), dof_handler, nodes);

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


        IPeriodicityMapper* pm = nullptr;
        if (per_type == "periodic-line" or per_type == "anti-periodic-line")
            pm = new LinePeriodicityMapper{dofs_1, dofs_2, dof_to_node};
        else if(per_type == "periodic-circle" or per_type == "anti-periodic-circle")
            pm = new CirclePeriodicityMapper{dofs_1, dofs_2, dof_to_node};

        auto matched_pairs = pm->get_matched_pair_indices();

        for (auto matched_pair: matched_pairs){

            auto first = matched_pair.first;
            auto second = matched_pair.second;

            constraints.add_line(first);
            constraints.add_entry(first, second, weight);
        }
    }

    constraints.print(std::cout);
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
    Point<dim> a{0, 0};
}

template<int dim>
void LinearSolver<dim>::assemble_system() {

    WorkStream::run(dof_handler.begin_active(),
                    dof_handler.end(),
                    *this,
                    &LinearSolver::local_assemble_system,
                    &LinearSolver::copy_local_to_global,
                    AssemblyScratchData(fe),
                    AssemblyCopyData());
}

template<int dim>
LinearSolver<dim>::AssemblyScratchData::AssemblyScratchData(const FiniteElement<dim> &fe):
    fe_values(fe, QGauss<dim>(fe.degree + 1), update_values | update_gradients | update_quadrature_points | update_JxW_values),
    rhs_values(fe_values.get_quadrature().size())
    {}

template<int dim>
LinearSolver<dim>::AssemblyScratchData::AssemblyScratchData(const AssemblyScratchData& scratch_data):
    fe_values(scratch_data.fe_values.get_fe(), scratch_data.fe_values.get_quadrature(), update_values | update_gradients | update_quadrature_points | update_JxW_values),
    rhs_values(scratch_data.rhs_values.size())
{}

template<int dim>
void LinearSolver<dim>::local_assemble_system(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                              LinearSolver::AssemblyScratchData &scratch_data,
                                              LinearSolver::AssemblyCopyData &copy_data) {

    double nu = 0;
    double f = 0;
    Tensor<1, dim> Hc({0, 0});

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points = scratch_data.fe_values.get_quadrature().size();

    copy_data.cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
    copy_data.cell_rhs.reinit(dofs_per_cell);
    copy_data.local_dof_indices.resize(dofs_per_cell);

    scratch_data.fe_values.reinit(cell);

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
                copy_data.cell_rhs(i) += f*                         // f at cell
                        sd.fe_values.shape_value(i, q)*            // phi_i(x_q)
                        sd.fe_values.JxW(q);                        // dx

            // Magnet contribution
                copy_data.cell_rhs(i) += Hc*sd.fe_values.shape_grad(i, q)*sd.fe_values.JxW(q);
        }

    cell->get_dof_indices(copy_data.local_dof_indices);
}

template<int dim>
void LinearSolver<dim>::copy_local_to_global(const AssemblyCopyData &copy_data) {
    constraints.distribute_local_to_global(
            copy_data.cell_matrix,
            copy_data.cell_rhs,
            copy_data.local_dof_indices,
            system_matrix,
            system_rhs);
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

