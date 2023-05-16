//
// Created by gordan on 5/16/23.
//

#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>

#include "MagneticFluxPostprocessor.h"

using namespace dealii;

template class MagneticFluxPostprocessor<2>;


template <int dim>
MagneticFluxPostprocessor<dim>::MagneticFluxPostprocessor(const Triangulation<dim> &triangulation,
                          const Vector<double>& solution, const FE_Q<dim>& fe) {
    triangulation_ptr = &triangulation;
    solution_ptr = &solution;
    fe_ptr = &fe;
}

template<int dim>
void MagneticFluxPostprocessor<dim>::process() {

    DoFHandler<dim> dof_handler(*triangulation_ptr);
    dof_handler.distribute_dofs(*fe_ptr);

    std::vector<double> Bx;
    std::vector<double> By;

    QGauss<dim> quadrature_formula(fe_ptr->degree +1);

    FEValues<dim> fe_values(*fe_ptr, quadrature_formula, update_values | update_gradients | update_quadrature_points |
                                                    update_JxW_values);
    std::vector<Tensor<1, dim>> solution_gradients(quadrature_formula.size());
    for (auto cell : dof_handler.active_cell_iterators()){
        fe_values.reinit(cell);
        fe_values.get_function_gradients(*solution_ptr, solution_gradients);

        Bx.push_back(solution_gradients[0][1]);
        By.push_back(-solution_gradients[0][0]);
    }


}

