//
// Created by gordan on 5/16/23.
//

#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>

#include "include/MagneticEnergyDensityPostprocessor.h"

using namespace dealii;

template class MagneticEnergyDensityPostprocessor<2>;


template <int dim>
MagneticEnergyDensityPostprocessor<dim>::MagneticEnergyDensityPostprocessor(const std::unordered_map<int, double>& nu_map) {
    this->nu_map_ptr = &nu_map;
}

template<int dim>
void MagneticEnergyDensityPostprocessor<dim>::process(const Triangulation<dim>&  triangulation,
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

        double nu_cell = (*nu_map_ptr).at(cell->material_id());

        double e_cell = 0;
        for (const unsigned int q: fe_values.quadrature_point_indices())
        {
            double b_x = solution_gradients[q][1];
            double b_y = -solution_gradients[q][0];

            e_cell += (b_x*b_x+b_y*b_y)/2.0 * nu_cell;

        }
        result.push_back(e_cell);
    }


}

