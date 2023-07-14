//
// Created by gordan on 5/16/23.
//
#include <variant>

#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>

#include "processors/FluxLinkageScalarPostprocessor.h"
#include "ConstFSource.h"

using namespace dealii;

template class FluxLinkageScalarPostprocessor<2>;

template <int dim>
FluxLinkageScalarPostprocessor<dim>::FluxLinkageScalarPostprocessor(const unsigned int mat_id, std::unordered_map<int, std::variant<FSource*, std::pair<double, double>>>& f_map) {
    this->f_map_ptr = &f_map;
    this->mat_id = mat_id;
}

template<int dim>
void FluxLinkageScalarPostprocessor<dim>::process(const Triangulation<dim>&  triangulation,
                                                 const Vector<double>&      solution,
                                                 const FE_Q<dim>&           fe,
                                                 double& result) {

    this->triangulation_ptr = &triangulation;
    this->solution_ptr = &solution;
    this->fe_ptr = &fe;

    std::vector<Point<dim>> q_points;
    DoFHandler<dim> dof_handler(*this->triangulation_ptr);
    dof_handler.distribute_dofs(*this->fe_ptr);
    QGauss<dim> quadrature_formula(this->fe_ptr->degree +1);
    FEValues<dim> fe_values(*this->fe_ptr, quadrature_formula, update_values | update_gradients | update_quadrature_points |
                                                         update_JxW_values);
    std::vector<Tensor<1, dim>> solution_gradients(quadrature_formula.size());
    std::vector<double> solution_at_cell(quadrature_formula.size());

    result = 0;
    FSource* f_source;
    ConstFSource f_zero{0};
    double i_current = 0;
    for (auto& cell : dof_handler.active_cell_iterators()){
        if (cell->material_id() == mat_id){
            if (f_map_ptr){
                auto f_variant = (*f_map_ptr).at(cell->material_id());
                if(std::holds_alternative<FSource*>(f_variant))
                    f_source = std::get<FSource*>(f_variant);
            }else{
                f_source = &f_zero;
            }

            fe_values.reinit(cell);
            fe_values.get_function_gradients(*this->solution_ptr, solution_gradients);
            fe_values.get_function_values(*this->solution_ptr, solution_at_cell);
            q_points = fe_values.get_quadrature_points();

            for (auto q: fe_values.quadrature_point_indices()){

                auto x = q_points[q][0];
                auto y = q_points[q][1];
                auto J = f_source->get_value(x, y);
                auto u = solution_at_cell[q];
                auto JxW = fe_values.JxW(q);
                result += u*J*JxW;
                i_current += J*JxW;
            }
        }
    }

    result = result/i_current;
}


