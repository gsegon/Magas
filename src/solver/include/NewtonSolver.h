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

#ifndef SOLVER_NEWTONSOLVER_H
#define SOLVER_NEWTONSOLVER_H

#include <any>
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
struct NuHistory
{
    //TODO: Generalize
    double nu[4];
};

template <int dim>
class NewtonSolver : public Solver<dim>{

public:

    void setup_system(const bool initial_step);
    void solve(const double alpha);
    double compute_residual(double) const;
    void run();
    void local_assemble_system(const typename DoFHandler<dim>::active_cell_iterator& cell,
                               typename Solver<dim>::AssemblyScratchData& scratch,
                               typename Solver<dim>::AssemblyCopyData& copy_data);

private:


    Vector<double> current_solution;
    Vector<double> newton_update;

    bool is_initial = true;

};

#endif //SOLVER_NEWTONSOLVER_H
