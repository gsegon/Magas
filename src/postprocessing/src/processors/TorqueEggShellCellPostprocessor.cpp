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

#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>

#include "processors/TorqueEggShellCellPostprocessor.h"

using namespace dealii;

template class TorqueEggShellCellPostprocessor<2>;


template <int dim>
TorqueEggShellCellPostprocessor<dim>::TorqueEggShellCellPostprocessor(unsigned int egg_id, unsigned int shell_id) {
    this->egg_id = egg_id;
    this->shell_id = shell_id;
}

template<int dim>
void TorqueEggShellCellPostprocessor<dim>::process(const Triangulation<dim>&  triangulation,
                                                   const Vector<double>&      solution,
                                                   const FE_Q<dim>&           fe,
                                                   std::vector<double>& result) {

    triangulation_ptr = &triangulation;
    solution_ptr = &solution;
    fe_ptr = &fe;

    DoFHandler<dim> dof_handler(*triangulation_ptr);
    dof_handler.distribute_dofs(*fe_ptr);
    QGauss<dim> quadrature_formula(fe_ptr->degree + 1);
    FEValues<dim> fe_values(*fe_ptr, quadrature_formula, update_values | update_gradients | update_quadrature_points |
                                                         update_JxW_values);
    std::vector<Tensor<1, dim>> solution_gradients(quadrature_formula.size());

    std::set<int> mask_vertex_indices;
    std::set<int> mask_cell_indices;
    for (auto cell: dof_handler.active_cell_iterators()) {
        if (cell->material_id() == egg_id){
            for (unsigned int i = 0; i < (cell->n_faces()); i++) {
                if (cell->neighbor_index(i) != -1) {
                    auto neighbor_cell = cell->neighbor(i);
                    if (neighbor_cell->material_id() == shell_id) {
                        mask_cell_indices.insert(neighbor_cell->index());
                        auto face = cell->face(i);
                        for (unsigned int j = 0; j < face->n_vertices(); j++) {
                            mask_vertex_indices.insert(face->vertex_index(j));
                        }
                    }
                }
            }
        }
    }

    for (auto cell: dof_handler.active_cell_iterators()) {
        if (cell->material_id() == shell_id){
            for (unsigned int j=0; j<cell->n_vertices(); j++){
                if (std::count(mask_vertex_indices.begin(), mask_vertex_indices.end(), cell->vertex_index(j))){
                    mask_cell_indices.insert(cell->index());
                }
            }
        }
    }


    double nu_0 = 795774.715025;
    for (auto cell: dof_handler.active_cell_iterators()) {
        double cell_torque = 0;
        if(std::count(mask_cell_indices.begin(), mask_cell_indices.end(), cell->index())){
            fe_values.reinit(cell);
            fe_values.get_function_gradients(solution, solution_gradients);
            auto q_points = fe_values.get_quadrature_points();
            for (unsigned int j=0; j<cell->n_vertices(); j++){
                if (std::count(mask_vertex_indices.begin(), mask_vertex_indices.end(), cell->vertex_index(j))){
                    for (auto q : fe_values.quadrature_point_indices()){
                        auto x = q_points[q][0];
                        auto y = q_points[q][1];
                        auto Bx = solution_gradients[q][1];
                        auto By = -solution_gradients[q][0];
                        auto JxW = fe_values.JxW(q);

                        Tensor<1, dim> r({x, y});
                        auto gamma = fe_values.shape_grad(j, q);
                        auto cross_1 = r[0]*gamma[1] -r[1]*gamma[0];
                        auto cross_2 = r[0]*By -r[1]*Bx;

                        double term1 = nu_0*0.5*solution_gradients[q].norm_square()*cross_1;
                        double term2 = -nu_0*(Bx*gamma[0] + By*gamma[1])*cross_2;
                        cell_torque += (term1 + term2)*JxW;
                    }
                }
            }
        }
        result.push_back(cell_torque);
    }
}

