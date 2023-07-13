//
// Created by gordan on 5/16/23.
//

#ifndef SOLVER_EXPRESSIONCELLPOSTPROCESSOR_H
#define SOLVER_EXPRESSIONCELLPOSTPROCESSOR_H

#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/grid_in.h>
#include "CellPostprocessor.h"
#include "NuCurve.h"
#include "FSource.h"

using namespace dealii;

typedef std::unordered_map<int, NuCurve*> t_nu_map;
typedef std::unordered_map<int, std::variant<FSource*, std::pair<double, double>>> t_f_map;
typedef std::unordered_map<int, double> t_dc_map;
typedef std::unordered_map<std::string, std::vector<unsigned int>> t_per_map;
typedef std::unordered_map<std::string, double> t_cli_source_map;

template <int dim>
class ExpressionCellPostprocessor : public CellPostprocessor<dim> {

public:
    ExpressionCellPostprocessor(const std::string&);
    ExpressionCellPostprocessor(const std::string&, t_nu_map&);
    ExpressionCellPostprocessor(const std::string&, t_f_map&);
    ExpressionCellPostprocessor(const std::string&, t_nu_map&, t_f_map&);

    void process(const Triangulation<dim>&  triangulation,
                 const Vector<double>&      solution,
                 const FE_Q<dim>&           fe,
                 std::vector<double>& result);

private:
    const Triangulation<dim> *triangulation_ptr = nullptr;
    const Vector<double> *solution_ptr = nullptr;
    const FE_Q<dim> *fe_ptr = nullptr;
    const t_nu_map* nu_map_ptr = nullptr;
    t_f_map* f_map_ptr = nullptr;

    std::string user_expression;

};

#endif //SOLVER_EXPRESSIONCELLPOSTPROCESSOR_H