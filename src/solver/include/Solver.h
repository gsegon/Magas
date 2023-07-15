//
// Created by gordan on 7/14/23.
//

#ifndef MAGAS_SOLVER_H
#define MAGAS_SOLVER_H



#include <string>
#include <unordered_map>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

#include "NuCurveFactory.h"
#include "FSourceFactory.h"

using namespace dealii;

typedef std::unordered_map<int, NuCurve*> t_nu_map;
typedef std::unordered_map<int, std::variant<FSource*, std::pair<double, double>>> t_f_map;
typedef std::unordered_map<int, double> t_dc_map;
typedef std::unordered_map<std::string, std::vector<unsigned int>> t_per_map;
typedef std::map<std::string, std::string> t_postprocessor_strings;

template<int dim>
class Solver{
public:
    virtual void read_mesh(const std::string&) = 0;
    virtual void set_nu_map(t_nu_map) = 0;
    virtual void set_f_map(t_f_map) = 0;
    virtual void set_dc_map(t_dc_map) = 0;
    virtual void set_per_map(t_per_map) = 0;
    virtual void run() = 0;

    virtual Triangulation<dim>& get_triangulation() = 0;
    virtual Vector<double>& get_solution() = 0;
    virtual Vector<double>& get_rhs() = 0;
    virtual FE_Q<dim>& get_fe() = 0;
};

#endif //MAGAS_SOLVER_H