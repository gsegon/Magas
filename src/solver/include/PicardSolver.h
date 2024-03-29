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

#ifndef SOLVER_PICARDSOLVER_H
#define SOLVER_PICARDSOLVER_H

#include <any>

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


using namespace dealii;

template <int dim>
struct NuHistory
{
    //TODO: Generalize
    double nu[4];
};

template <int dim>
class PicardSolver{

public:
    PicardSolver();
    void read_mesh(const std::string&);
    void setup_system();
    void reinit_system();
    void assemble_system();
    void solve();
    void solve_nonlinear(int);
    void set_nu_map(std::unordered_map<int, std::any>);
    void set_f_map(std::unordered_map<int, double>);
    void set_dc_map(std::unordered_map<int, double>);
    void setup_cell_nu_history();
    void update_cell_nu_history();
    void initialize_cell_nu_history(double);

    Triangulation<dim>& get_triangulation();
    Vector<double>& get_solution();
    Vector<double>& get_rhs();
    FE_Q<dim>& get_fe();

private:

    Triangulation<dim> triangulation;
    FE_Q<dim> fe;
    DoFHandler<dim> dof_handler;
    QGauss<dim> quadrature_formula;

    SparsityPattern sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double> solution;
    Vector<double> system_rhs;

    std::unordered_map<int, std::any> nu_map;
    std::unordered_map<int, double> f_map;
    std::unordered_map<int, double> dc_map;
    std::vector<NuHistory<dim>> quadrature_nu_history;

};




#endif //SOLVER_PICARDSOLVER_H
