//
// Created by gordan on 5/16/23.
//

#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>

#include "MatIDPostprocessor.h"

using namespace dealii;

template class MatIDPostprocessor<2>;


template <int dim>
MatIDPostprocessor<dim>::MatIDPostprocessor() {
}

template<int dim>
void MatIDPostprocessor<dim>::process(const Triangulation<dim>&  triangulation,
                                             const Vector<double>&      solution,
                                             const FE_Q<dim>&           fe,
                                             std::vector<double>& result) {

    triangulation_ptr = &triangulation;
    solution_ptr = &solution;
    fe_ptr = &fe;

    DoFHandler<dim> dof_handler(*triangulation_ptr);
    dof_handler.distribute_dofs(*fe_ptr);

    auto mat_id = &result;

    QGauss<dim> quadrature_formula(fe_ptr->degree +1);

    FEValues<dim> fe_values(*fe_ptr, quadrature_formula, update_values | update_gradients | update_quadrature_points |
                                                    update_JxW_values);
    std::vector<Tensor<1, dim>> solution_gradients(quadrature_formula.size());
    for (auto cell : dof_handler.active_cell_iterators()){
       mat_id->push_back(cell->material_id());
    }



}

