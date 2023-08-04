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
ArkkioScalarPostprocessor<dim>::ArkkioScalarPostprocessor(const unsigned int mat_id, const std::unordered_map<int, NuCurve*>& nu_map) {
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
    double nu_q1 = 0;
    double nu_q2 = 0;
    double nu_q3 = 0;
    double nu_q4 = 0;
    double nu_0 = 795774.715025;
    if (nu_map_ptr){
        NuCurve* bh = ((*nu_map_ptr).at(mat_id));
        nu_q1 = bh->get_nu(solution_gradients[0].norm());
        nu_q2 = bh->get_nu(solution_gradients[1].norm());
        nu_q3 = bh->get_nu(solution_gradients[2].norm());
        nu_q4 = bh->get_nu(solution_gradients[3].norm());
    }

    AssertThrow(std::abs(nu_q1-nu_0) < 1e-3, ExcInternalError())
    AssertThrow(std::abs(nu_q2-nu_0) < 1e-3, ExcInternalError())
    AssertThrow(std::abs(nu_q3-nu_0) < 1e-3, ExcInternalError())
    AssertThrow(std::abs(nu_q4-nu_0) < 1e-3, ExcInternalError())

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

    result *= nu_q1/(r2-r1);

}


