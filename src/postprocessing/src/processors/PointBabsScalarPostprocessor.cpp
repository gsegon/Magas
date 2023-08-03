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
#include <deal.II/numerics/vector_tools_point_value.h>
#include <deal.II/numerics/vector_tools_point_gradient.h>
#include "exprtk.hpp"

#include "processors/PointBabsScalarPostprocessor.h"

typedef exprtk::symbol_table<double> symbol_table_t;
typedef exprtk::expression<double>   expression_t;
typedef exprtk::parser<double>       parser_t;

using namespace dealii;

template class PointBabsScalarPostprocessor<2>;



template <int dim>
PointBabsScalarPostprocessor<dim>::PointBabsScalarPostprocessor(std::vector<double> point) {
    this->point = point;
}


template<int dim>
void PointBabsScalarPostprocessor<dim>::process(const Triangulation<dim>&  triangulation,
                                                 const Vector<double>&      solution,
                                                 const FE_Q<dim>&           fe,
                                                 double& result) {

    this->triangulation_ptr = &triangulation;
    this->solution_ptr = &solution;
    this->fe_ptr = &fe;

    std::vector<Point<dim>> q_points;
    DoFHandler<dim> dof_handler(*(this->triangulation_ptr));
    dof_handler.distribute_dofs(*(this->fe_ptr));

    Point<dim> deal_point{point[0], point[1]};

    result = 0;
    auto value_gradient = VectorTools::point_gradient(dof_handler, solution, deal_point);
    result = value_gradient.norm();

}


