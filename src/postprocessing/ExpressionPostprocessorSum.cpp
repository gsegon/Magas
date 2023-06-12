//
// Created by gordan on 5/16/23.
//

#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include "exprtk.hpp"

#include "include/ExpressionPostprocessorSum.h"

typedef exprtk::symbol_table<double> symbol_table_t;
typedef exprtk::expression<double>   expression_t;
typedef exprtk::parser<double>       parser_t;

using namespace dealii;

template class ExpressionPostprocessorSum<2>;


template <int dim>
ExpressionPostprocessorSum<dim>::ExpressionPostprocessorSum(const std::string& user_expr, const std::unordered_map<int, double>& nu_map) {
    user_expression = user_expr;
    nu_map_ptr = &nu_map;
}


template<int dim>
void ExpressionPostprocessorSum<dim>::process(const Triangulation<dim>&  triangulation,
                                             const Vector<double>&      solution,
                                             const FE_Q<dim>&           fe,
                                             double& result) {

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

    double JxW_q1 = 0;
    double JxW_q2 = 0;
    double JxW_q3 = 0;
    double JxW_q4 = 0;

    double nu = 0;

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

    symbol_table.add_variable("mat_id", mat_id);
    symbol_table.add_variable("nu", nu);

    symbol_table.add_constants();
    expression.register_symbol_table(symbol_table);
    parser.compile(user_expression, expression);


    DoFHandler<dim> dof_handler(*triangulation_ptr);
    dof_handler.distribute_dofs(*fe_ptr);

    QGauss<dim> quadrature_formula(fe_ptr->degree +1);
    FEValues<dim> fe_values(*fe_ptr, quadrature_formula, update_values | update_gradients | update_quadrature_points |
                                                         update_JxW_values);
    std::vector<Tensor<1, dim>> solution_gradients(quadrature_formula.size());

    result = 0;
    for (auto& cell : dof_handler.active_cell_iterators()){

        fe_values.reinit(cell);
        fe_values.get_function_gradients(*solution_ptr, solution_gradients);

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

        nu = (*nu_map_ptr).at(cell->material_id());

        result += expression.value();

    }

}


