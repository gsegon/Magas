//
// Created by gordan on 5/16/23.
//

#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>

#include "processors/EggShellCellPostprocessor.h"

using namespace dealii;

template class EggShellCellPostprocessor<2>;


template <int dim>
EggShellCellPostprocessor<dim>::EggShellCellPostprocessor(unsigned int egg_id, unsigned int shell_id) {
    this->egg_id = egg_id;
    this->shell_id = shell_id;
}

template<int dim>
void EggShellCellPostprocessor<dim>::process(const Triangulation<dim>&  triangulation,
                                             const Vector<double>&      solution,
                                             const FE_Q<dim>&           fe,
                                             std::vector<double>& result) {

    triangulation_ptr = &triangulation;
    solution_ptr = &solution;
    fe_ptr = &fe;

    DoFHandler<dim> dof_handler(*triangulation_ptr);
    dof_handler.distribute_dofs(*fe_ptr);

    std::vector<double> *mask = &result;

    QGauss<dim> quadrature_formula(fe_ptr->degree + 1);

    FEValues<dim> fe_values(*fe_ptr, quadrature_formula, update_values | update_gradients | update_quadrature_points |
                                                         update_JxW_values);
    std::vector<Tensor<1, dim>> solution_gradients(quadrature_formula.size());
    std::set<int> mask_vertex_indices;
    std::set<int> mask_cell_indices;
    double different = 0;
    for (auto cell: dof_handler.active_cell_iterators()) {
        if (cell->material_id() == egg_id){
            for (unsigned int i = 0; i < (cell->n_faces()); i++) {
                if (cell->neighbor_index(i) != -1) {
                    auto neighbor_cell = cell->neighbor(i);
                    if (neighbor_cell->material_id() == shell_id) {
                        mask_cell_indices.insert(neighbor_cell->index());
                        auto face = cell->face(i);
                        for (unsigned int j = 0; j < face->n_vertices(); j++) {
                            auto vertex = face->vertex(j);
                            mask_vertex_indices.insert(face->vertex_index(j));
                        }
                    }
                }
            }
        }
    }

    for (auto cell: dof_handler.active_cell_iterators()) {
        if (cell->material_id() == shell_id){
            for (int j=0; j<cell->n_vertices(); j++){
                if (std::count(mask_vertex_indices.begin(), mask_vertex_indices.end(), cell->vertex_index(j))){
                    mask_cell_indices.insert(cell->index());
                }
            }
        }
    }

    for (auto cell: dof_handler.active_cell_iterators())
        result.push_back(0);

    for (auto cell_index: mask_cell_indices) {
        result.insert(result.begin() + cell_index, 1);
    }
}

