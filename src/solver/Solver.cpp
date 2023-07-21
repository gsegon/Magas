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

#include "PeriodicityMapperFactory.h"
#include "Solver.h"
#include "ConstFSource.h"

//template class Solver<2>;

using namespace dealii;

template class Solver<2>;

template <int dim>
Solver<dim>::Solver(): fe(1), dof_handler(triangulation), quadrature_formula(fe.degree+1)  {}

template<int dim>
void Solver<dim>::read_mesh(const std::string& mesh_filepath) {
    GridIn<dim> grid_in;
    grid_in.attach_triangulation(triangulation);
    std::ifstream input_file(mesh_filepath);
    grid_in.read_msh(input_file);
}

template<int dim>
Triangulation<dim>& Solver<dim>::get_triangulation(){
    return this->triangulation;
}

template<int dim>
void Solver<dim>::assemble_system() {

    // Rest system matrix, system rhs before assembly
    system_matrix.reinit(sparsity_pattern);
    system_rhs.reinit(dof_handler.n_dofs());

    WorkStream::run(dof_handler.begin_active(),
                    dof_handler.end(),
                    *this,
                    &Solver::local_assemble_system,
                    &Solver::copy_local_to_global,
                    AssemblyScratchData(fe),
                    AssemblyCopyData());
}

template<int dim>
Solver<dim>::AssemblyScratchData::AssemblyScratchData(const FiniteElement<dim> &fe):
    fe_values(fe, QGauss<dim>(fe.degree + 1), update_values | update_gradients | update_quadrature_points | update_JxW_values),
    rhs_values(fe_values.get_quadrature().size())
    {}

template<int dim>
Solver<dim>::AssemblyScratchData::AssemblyScratchData(const AssemblyScratchData& scratch_data):
    fe_values(scratch_data.fe_values.get_fe(), scratch_data.fe_values.get_quadrature(), update_values | update_gradients | update_quadrature_points | update_JxW_values),
    rhs_values(scratch_data.rhs_values.size())
{}

template<int dim>
void Solver<dim>::copy_local_to_global(const AssemblyCopyData &copy_data) {
    constraints.distribute_local_to_global(
            copy_data.cell_matrix,
            copy_data.cell_rhs,
            copy_data.local_dof_indices,
            system_matrix,
            system_rhs);
}

template<int dim>
void Solver<dim>::set_nu_map(std::unordered_map<int, NuCurve*> map) {
    this->nu_map = map;
}

template<int dim>
void Solver<dim>::set_f_map(std::unordered_map<int, std::variant<FSource*, std::pair<double, double>>> map) {
    this->f_map = map;
}

template<int dim>
void Solver<dim>::set_dc_map(std::unordered_map<int, double> map) {
    this->dc_map = map;
}

template<int dim>
void Solver<dim>::set_per_map(std::unordered_map<std::string, std::vector<unsigned int>> map) {
    this->per_map = map;
}


template<int dim>
Vector<double>& Solver<dim>::get_solution(){
    return this->solution;
}

template<int dim>
Vector<double>& Solver<dim>::get_rhs(){
    return this->system_rhs;
}

template<int dim>
FE_Q<dim>& Solver<dim>::get_fe(){
    return this->fe;
}


