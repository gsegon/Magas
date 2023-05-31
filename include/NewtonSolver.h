//
// Created by gordan on 5/10/23.
//

#ifndef SOLVER_NEWTONSOLVER_H
#define SOLVER_NEWTONSOLVER_H

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
class NewtonSolver{

public:
    NewtonSolver();
    void read_mesh(const std::string&);
    void setup_system(const bool initial_step);
    void assemble_system();
    void solve(const double alpha);
//    void solve_nonlinear(int);
    void set_nu_map(std::unordered_map<int, std::any>);
    void set_f_map(std::unordered_map<int, std::variant<double, std::pair<double, double>>>);
    void set_dc_map(std::unordered_map<int, double>);
    double compute_residual() const;

    Triangulation<dim>& get_triangulation();
    Vector<double>& get_solution();
    Vector<double>& get_current_solution();
    Vector<double>& get_rhs();
    FE_Q<dim>& get_fe();

private:

    Triangulation<dim> triangulation;
    FE_Q<dim> fe;
    DoFHandler<dim> dof_handler;
    QGauss<dim> quadrature_formula;

    SparsityPattern sparsity_pattern;
    SparseMatrix<double> system_matrix;

    AffineConstraints<double> hanging_node_constraints;

    Vector<double> solution;
    Vector<double> system_rhs;

    Vector<double> current_solution;
    Vector<double> newton_update;

    std::unordered_map<int, std::any> nu_map;
    std::unordered_map<int, std::variant<double, std::pair<double, double>>> f_map;
    std::unordered_map<int, double> dc_map;

};

#endif //SOLVER_NEWTONSOLVER_H
