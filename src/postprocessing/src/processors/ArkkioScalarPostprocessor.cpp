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
#include "exprtk.hpp"

#include "processors/ArkkioScalarPostprocessor.h"

typedef exprtk::symbol_table<double> symbol_table_t;
typedef exprtk::expression<double>   expression_t;
typedef exprtk::parser<double>       parser_t;

using namespace dealii;

template class ArkkioScalarPostprocessor<2>;



template <int dim>
ArkkioScalarPostprocessor<dim>::ArkkioScalarPostprocessor(const unsigned int mat_id, const std::unordered_map<int, double>& nu_map) {
    this->nu_map_ptr = &nu_map;
    this->mat_id = mat_id;
}


template<int dim>
void ArkkioScalarPostprocessor<dim>::process(const Triangulation<dim>&  triangulation,
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

    // Checks
    double nu = 0;
    double nu_0 = 795774.715025;
    if (nu_map_ptr)
        nu = (*nu_map_ptr).at(mat_id);

    AssertThrow(std::abs(nu-nu_0) < 1e-3, ExcInternalError())

    double r1=std::numeric_limits<double>::infinity();
    double r2=-std::numeric_limits<double>::infinity();
    result = 0;
    for (auto& cell : dof_handler.active_cell_iterators()){
        if (cell->material_id() == mat_id){
            fe_values.reinit(cell);
            fe_values.get_function_gradients(*this->solution_ptr, solution_gradients);
            fe_values.get_function_values(*this->solution_ptr, solution_at_cell);
            q_points = fe_values.get_quadrature_points();

            for (auto vertex_index : GeometryInfo<dim>::vertex_indices()){
                double r = cell->vertex(vertex_index).norm();
                if(r < r1) r1 = r;
                if(r > r2) r2 = r;
            }
            for (auto q: fe_values.quadrature_point_indices()){
                auto x = q_points[q][0];
                auto y = q_points[q][1];
                auto Bx = solution_gradients[q][1];
                auto By = -solution_gradients[q][0];
                auto JxW = fe_values.JxW(q);
                result += (Bx*By*x*x +(By*By-Bx*Bx)*x*y-Bx*By*y*y)/sqrt(x*x+y*y)*JxW;
            }
        }
    }

    result *= nu/(r2-r1);

}


