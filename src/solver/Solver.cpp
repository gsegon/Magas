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

using namespace dealii;

template<int dim>
void set_rotating_bound_data(unsigned int domain1,
                             unsigned int domain2,
                             const Triangulation<dim>&  triangulation,
                             const FE_Q<dim>&           fe,
                             std::vector<unsigned int>& cell_indices,
                             std::vector<unsigned int>& dofs
){


    DoFHandler<dim> dof_handler(triangulation);
    dof_handler.distribute_dofs(fe);

    QGauss<dim> quadrature_formula(fe.degree + 1);
    FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_gradients | update_quadrature_points |
                                                    update_JxW_values);
    std::vector<Tensor<1, dim>> solution_gradients(quadrature_formula.size());

    // Stator
    std::set<int> mask_stator_vertex_indices;
    std::set<int> mask_stator_cell_indices;
    std::set<std::pair<unsigned int, int>> stator_dof_to_vertex_index;
    std::set<unsigned int> dofs1;

    for (auto cell: dof_handler.active_cell_iterators()) {
        if (cell->material_id() == domain1){
            for (unsigned int i = 0; i < (cell->n_faces()); i++) {
                if (cell->neighbor_index(i) != -1) {
                    auto neighbor_cell = cell->neighbor(i);
                    if (neighbor_cell->material_id() == domain2) {
                        mask_stator_cell_indices.insert(neighbor_cell->index());
                        auto face = cell->face(i);
                        std::vector< types::global_dof_index > global_dof_indices(2);
                        face->get_dof_indices(global_dof_indices);
                        for (unsigned int j = 0; j < face->n_vertices(); j++) {
                            auto vertex = face->vertex(j);
                            mask_stator_vertex_indices.insert(face->vertex_index(j));
                            stator_dof_to_vertex_index.insert({global_dof_indices[j], face->vertex_index(j)});
                            dofs1.insert(global_dof_indices[j]);
                        }
                    }
                }
            }
        }
    }

    for (auto cell: dof_handler.active_cell_iterators())
        if (cell->material_id() == domain2)
            for (unsigned int j=0; j<cell->n_vertices(); j++)
                if (std::count(mask_stator_vertex_indices.begin(), mask_stator_vertex_indices.end(), cell->vertex_index(j)))
                    mask_stator_cell_indices.insert(cell->index());

    for (auto index: mask_stator_cell_indices)
        cell_indices.push_back(index);

    for (auto dof: dofs1)
        dofs.push_back(dof);

}


template class Solver<2>;

template <int dim>
Solver<dim>::Solver(): fe(1), dof_handler(triangulation), quadrature_formula(fe.degree+1)  {}


template<int dim>
void Solver<dim>::setup_rotation(unsigned int a, unsigned int b, int offset) {

//    dof_handler.distribute_dofs(fe);
    set_rotating_bound_data(a, b, triangulation, fe, rot_cell_indices, rot_dofs);

    std::vector<Point<dim>> nodes(dof_handler.n_dofs());
    DoFTools::map_dofs_to_support_points(MappingQ1<dim>(), dof_handler, nodes);
    std::map<unsigned int, std::vector<double>> dof_to_node;
    for (auto dof : rot_dofs){
        dof_to_node[dof] = {nodes[dof][0], nodes[dof][1]};
    }
    sr = new SlidingRotation{rot_dofs, dof_to_node, offset};
    std::cout << "Sliding boundary has " << sr->get_dofs().size() << " dofs." << std::endl;

}

template<int dim>
void Solver<dim>::extend_dsp(DynamicSparsityPattern& dsp){
    std::vector<types::global_dof_index> local_dof_indices(4);
    for (auto cell: dof_handler.active_cell_iterators()){
        cell->get_dof_indices(local_dof_indices);

        if (std::count(rot_cell_indices.begin(), rot_cell_indices.end(), cell->index())){
            for (auto& local_dof_index : local_dof_indices){
                local_dof_index = sr->get_mapped(local_dof_index);
            }
            for (auto i : local_dof_indices)
                for (auto j : local_dof_indices){
                    dsp.add(i, j);
                }
        }
    }
}

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
void Solver<dim>::set_nu_map(t_nu_map map) {
    this->nu_map = map;
}

template<int dim>
void Solver<dim>::set_f_map(t_f_map map) {
    this->f_map = map;
}

template<int dim>
void Solver<dim>::set_dc_map(t_dc_map map) {
    this->dc_map = map;
}

template<int dim>
void Solver<dim>::set_per_map(t_per_map map) {
    this->per_map = map;
}

template<int dim>
void Solver<dim>::set_rot_map(t_rot_map map) {
    this->rot_map = map;
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