//
// Created by gordan on 5/10/23.
//

#ifndef SOLVER_LINEARSOLVER_H
#define SOLVER_LINEARSOLVER_H

#include <variant>

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
#include <deal.II/numerics/data_out.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/manifold_lib.h>
#include "NuCurve.h"
#include "FSource.h"
#include "Solver.h"

using namespace dealii;

template <int dim>
class LinearSolver: public Solver<dim>{

public:
//    LinearSolver();
    void setup_system();
    void solve();
    void run();
    void local_assemble_system(const typename DoFHandler<dim>::active_cell_iterator& cell,
                               typename Solver<dim>::AssemblyScratchData& scratch,
                               typename Solver<dim>::AssemblyCopyData& copy_data);

};

#endif //SOLVER_LINEARSOLVER_H
