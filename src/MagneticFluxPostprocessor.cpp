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
MagneticFluxPostprocessor<dim>::MagneticFluxPostprocessor() {
    this->abs = true;
}


template <int dim>
MagneticFluxPostprocessor<dim>::MagneticFluxPostprocessor(unsigned int component) {
    this->component = component;
}

template<int dim>
void MagneticFluxPostprocessor<dim>::process(const Triangulation<dim>&  triangulation,
                                             const Vector<double>&      solution,
                                             const FE_Q<dim>&           fe,
                                             std::vector<double>& result) {

    triangulation_ptr = &triangulation;
    solution_ptr = &solution;
    fe_ptr = &fe;

    DoFHandler<dim> dof_handler(*triangulation_ptr);
    dof_handler.distribute_dofs(*fe_ptr);

    QGauss<dim> quadrature_formula(fe_ptr->degree +1);

    FEValues<dim> fe_values(*fe_ptr, quadrature_formula, update_values | update_gradients | update_quadrature_points |
                                                    update_JxW_values);
    std::vector<Tensor<1, dim>> solution_gradients(quadrature_formula.size());
    for (auto cell : dof_handler.active_cell_iterators()){
        fe_values.reinit(cell);
        fe_values.get_function_gradients(*solution_ptr, solution_gradients);

        if (this->abs) {
            result.push_back(std::sqrt(std::pow(solution_gradients[0][1], 2) + std::pow(-solution_gradients[0][0], 2)));
        }
        else{
            if (component == 0){
                result.push_back(solution_gradients[0][1]);
            }
            else if (component == 1){
                result.push_back(-solution_gradients[0][0]);
            }
        }



    }



}

