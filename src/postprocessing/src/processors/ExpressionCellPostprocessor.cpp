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

#include "processors/ExpressionCellPostprocessor.h"
#include "NuCurve.h"
#include "FSource.h"
#include "ConstFSource.h"

typedef exprtk::symbol_table<double> symbol_table_t;
typedef exprtk::expression<double>   expression_t;
typedef exprtk::parser<double>       parser_t;

using namespace dealii;

template class ExpressionCellPostprocessor<2>;


template <int dim>
ExpressionCellPostprocessor<dim>::ExpressionCellPostprocessor(const std::string& user_expr) {
    user_expression = user_expr;
}

template <int dim>
ExpressionCellPostprocessor<dim>::ExpressionCellPostprocessor(const std::string& user_expr, t_nu_map& nu_map) {
    user_expression = user_expr;
    nu_map_ptr = &nu_map;
}

template <int dim>
ExpressionCellPostprocessor<dim>::ExpressionCellPostprocessor(const std::string& user_expr, t_f_map& f_map) {
    user_expression = user_expr;
    f_map_ptr = &f_map;
}

template <int dim>
ExpressionCellPostprocessor<dim>::ExpressionCellPostprocessor(const std::string& user_expr, t_nu_map& nu_map, t_f_map& f_map) {
    user_expression = user_expr;
    nu_map_ptr = &nu_map;
    f_map_ptr = &f_map;
}


template<int dim>
void ExpressionCellPostprocessor<dim>::process(const Triangulation<dim>&  triangulation,
                                               const Vector<double>&      solution,
                                               const FE_Q<dim>&           fe,
                                               std::vector<double>& result) {

    triangulation_ptr = &triangulation;
    solution_ptr = &solution;
    fe_ptr = &fe;

    symbol_table_t symbol_table;
    expression_t expression;
    parser_t parser;

    std::vector<Point<2>> q_points;
    double mat_id = -1;

    double x_q1 = 0;
    double x_q2 = 0;
    double x_q3 = 0;
    double x_q4 = 0;

    double y_q1 = 0;
    double y_q2 = 0;
    double y_q3 = 0;
    double y_q4 = 0;

    double Bx_q1 = 0;
    double Bx_q2 = 0;
    double Bx_q3 = 0;
    double Bx_q4 = 0;

    double By_q1 = 0;
    double By_q2 = 0;
    double By_q3 = 0;
    double By_q4 = 0;

    double u_q1 = 0;
    double u_q2 = 0;
    double u_q3 = 0;
    double u_q4 = 0;

    double JxW_q1 = 0;
    double JxW_q2 = 0;
    double JxW_q3 = 0;
    double JxW_q4 = 0;

    double J = 0;
    double J_q1 = 0;
    double J_q2 = 0;
    double J_q3 = 0;
    double J_q4 = 0;

    double nu_q1 = 0;
    double nu_q2 = 0;
    double nu_q3 = 0;
    double nu_q4 = 0;

    double e_q1 = 0;
    double e_q2 = 0;
    double e_q3 = 0;
    double e_q4 = 0;

    symbol_table.add_variable("mat_id", mat_id);

    symbol_table.add_variable("x_q1", x_q1);
    symbol_table.add_variable("x_q2", x_q2);
    symbol_table.add_variable("x_q3", x_q3);
    symbol_table.add_variable("x_q4", x_q4);

    symbol_table.add_variable("y_q1", y_q1);
    symbol_table.add_variable("y_q2", y_q2);
    symbol_table.add_variable("y_q3", y_q3);
    symbol_table.add_variable("y_q4", y_q4);

    symbol_table.add_variable("Bx_q1", Bx_q1);
    symbol_table.add_variable("Bx_q2", Bx_q2);
    symbol_table.add_variable("Bx_q3", Bx_q3);
    symbol_table.add_variable("Bx_q4", Bx_q4);

    symbol_table.add_variable("By_q1", By_q1);
    symbol_table.add_variable("By_q2", By_q2);
    symbol_table.add_variable("By_q3", By_q3);
    symbol_table.add_variable("By_q4", By_q4);

    symbol_table.add_variable("JxW_q1", JxW_q1);
    symbol_table.add_variable("JxW_q2", JxW_q2);
    symbol_table.add_variable("JxW_q3", JxW_q3);
    symbol_table.add_variable("JxW_q4", JxW_q4);

    symbol_table.add_variable("u_q1", u_q1);
    symbol_table.add_variable("u_q2", u_q2);
    symbol_table.add_variable("u_q3", u_q3);
    symbol_table.add_variable("u_q4", u_q4);

    if (f_map_ptr){
        symbol_table.add_variable("J_q1", J_q1);
        symbol_table.add_variable("J_q2", J_q2);
        symbol_table.add_variable("J_q3", J_q3);
        symbol_table.add_variable("J_q4", J_q4);
    }

    if (nu_map_ptr){
        symbol_table.add_variable("nu_q1", nu_q1);
        symbol_table.add_variable("nu_q2", nu_q2);
        symbol_table.add_variable("nu_q3", nu_q3);
        symbol_table.add_variable("nu_q4", nu_q4);

        symbol_table.add_variable("e_q1", e_q1);
        symbol_table.add_variable("e_q2", e_q2);
        symbol_table.add_variable("e_q3", e_q3);
        symbol_table.add_variable("e_q4", e_q4);
    }

    symbol_table.add_constants();
    expression.register_symbol_table(symbol_table);
    parser.compile(user_expression, expression);

    DoFHandler<dim> dof_handler(*triangulation_ptr);
    dof_handler.distribute_dofs(*fe_ptr);

    QGauss<dim> quadrature_formula(fe_ptr->degree +1);

    FEValues<dim> fe_values(*fe_ptr, quadrature_formula, update_values | update_gradients | update_quadrature_points |
                                                         update_JxW_values);
    std::vector<Tensor<1, dim>> solution_gradients(quadrature_formula.size());
    std::vector<double> solution_at_cell(quadrature_formula.size());
    FSource* f_source;
    ConstFSource f_zero{0};
    for (auto& cell : dof_handler.active_cell_iterators()){

        fe_values.reinit(cell);
        fe_values.get_function_gradients(*solution_ptr, solution_gradients);
        fe_values.get_function_values(*solution_ptr, solution_at_cell);

        q_points = fe_values.get_quadrature_points();

        mat_id = cell->material_id();
        Bx_q1 = solution_gradients[0][1];
        Bx_q2 = solution_gradients[1][1];
        Bx_q3 = solution_gradients[2][1];
        Bx_q4 = solution_gradients[3][1];

        By_q1 = -solution_gradients[0][0];
        By_q2 = -solution_gradients[1][0];
        By_q3 = -solution_gradients[2][0];
        By_q4 = -solution_gradients[3][0];

        x_q1 = q_points[0][0];
        x_q2 = q_points[1][0];
        x_q3 = q_points[2][0];
        x_q4 = q_points[3][0];

        y_q1 = q_points[0][1];
        y_q2 = q_points[1][1];
        y_q3 = q_points[2][1];
        y_q4 = q_points[3][1];

        JxW_q1 = fe_values.JxW(0);
        JxW_q2 = fe_values.JxW(1);
        JxW_q3 = fe_values.JxW(2);
        JxW_q4 = fe_values.JxW(3);

        u_q1 = solution_at_cell[0];
        u_q2 = solution_at_cell[1];
        u_q3 = solution_at_cell[2];
        u_q4 = solution_at_cell[3];

        if (f_map_ptr){
            auto f_variant = (*f_map_ptr).at(cell->material_id());
            if(std::holds_alternative<FSource*>(f_variant)){
                f_source = std::get<FSource*>(f_variant);
            }
            else
                f_source = &f_zero;

            J_q1 = f_source->get_value(x_q1, y_q1);
            J_q2 = f_source->get_value(x_q2, y_q2);
            J_q3 = f_source->get_value(x_q3, y_q3);
            J_q4 = f_source->get_value(x_q4, y_q4);
        }

        if (nu_map_ptr){
            NuCurve* bh = ((*nu_map_ptr).at(cell->material_id()));
            nu_q1 = bh->get_nu(solution_gradients[0].norm());
            nu_q2 = bh->get_nu(solution_gradients[1].norm());
            nu_q3 = bh->get_nu(solution_gradients[2].norm());
            nu_q4 = bh->get_nu(solution_gradients[3].norm());

            e_q1 = bh->get_energy(solution_gradients[0].norm());
            e_q2 = bh->get_energy(solution_gradients[1].norm());
            e_q3 = bh->get_energy(solution_gradients[2].norm());
            e_q4 = bh->get_energy(solution_gradients[3].norm());
        }
        result.push_back(expression.value());
    }
}


