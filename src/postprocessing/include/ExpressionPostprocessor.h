//
// Created by gordan on 5/16/23.
//

#ifndef SOLVER_EXPRESSIONPOSTPROCESSOR_H
#define SOLVER_EXPRESSIONPOSTPROCESSOR_H

#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/grid_in.h>
#include "Postprocessor.h"

using namespace dealii;

template <int dim>
class ExpressionPostprocessor : public Postprocessor<dim> {

public:
    ExpressionPostprocessor(const std::string&);
    ExpressionPostprocessor(const std::string&, const std::unordered_map<int, double>&);
    ExpressionPostprocessor(const std::string&, std::unordered_map<int, std::variant<double, std::pair<double, double>>>&);
    ExpressionPostprocessor(const std::string&, const std::unordered_map<int, double>&, std::unordered_map<int, std::variant<double, std::pair<double, double>>>&);

    void process(const Triangulation<dim>&  triangulation,
                 const Vector<double>&      solution,
                 const FE_Q<dim>&           fe,
                 std::vector<double>& result);

private:
    const Triangulation<dim> *triangulation_ptr = nullptr;
    const Vector<double> *solution_ptr = nullptr;
    const FE_Q<dim> *fe_ptr = nullptr;
    const std::unordered_map<int, double>* nu_map_ptr = nullptr;
    std::unordered_map<int, std::variant<double, std::pair<double, double>>>* f_map_ptr = nullptr;

    std::string user_expression;

};

#endif //SOLVER_EXPRESSIONPOSTPROCESSOR_H