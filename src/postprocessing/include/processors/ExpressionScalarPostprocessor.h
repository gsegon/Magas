//
// Created by gordan on 5/16/23.
//

#ifndef SOLVER_EXPRESSIONSCALARPOSTPROCESSOR_H
#define SOLVER_EXPRESSIONSCALARPOSTPROCESSOR_H

#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/grid_in.h>
#include "ScalarPostprocessor.h"

using namespace dealii;

template <int dim>
class ExpressionScalarPostprocessor : public ScalarPostprocessor<dim> {

public:
    explicit ExpressionScalarPostprocessor(const std::string&);
    ExpressionScalarPostprocessor(const std::string&, const std::unordered_map<int, double>&);
    ExpressionScalarPostprocessor(const std::string&, std::unordered_map<int, std::variant<double, std::pair<double, double>>>&);
    ExpressionScalarPostprocessor(const std::string&, const std::unordered_map<int, double>&, std::unordered_map<int, std::variant<double, std::pair<double, double>>>&);

    void process(const Triangulation<dim>&  triangulation,
                 const Vector<double>&      solution,
                 const FE_Q<dim>&           fe,
                 double& result);

private:

    std::string user_expression;
    const std::unordered_map<int, double>* nu_map_ptr = nullptr;
    std::unordered_map<int, std::variant<double, std::pair<double, double>>>* f_map_ptr = nullptr;

};

#endif //SOLVER_EXPRESSIONSCALARPOSTPROCESSOR_H