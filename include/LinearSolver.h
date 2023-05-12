//
// Created by gordan on 5/10/23.
//

#ifndef SOLVER_LINEARSOLVER_H
#define SOLVER_LINEARSOLVER_H

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
class LinearSolver{

public:
    LinearSolver();
    void read_mesh(std::string);
    void setup_system();
    void assemble_system();
    void solve();
    void output_results(const std::string filename) const;
    void set_nu_map(std::unordered_map<int, double> map);
    void set_f_map(std::unordered_map<int, double> map);
    void set_dc_map(std::unordered_map<int, double> map);
    Triangulation<dim>& get_triangulation();
//    void run();

private:

    Triangulation<dim> triangulation;
    FE_Q<dim> fe;
    DoFHandler<dim> dof_handler;
    QGauss<dim> quadrature_formula;

    SparsityPattern sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double> solution;
    Vector<double> system_rhs;
    std::unordered_map<int, double> nu_map;
    std::unordered_map<int, double> f_map;
    std::unordered_map<int, double> dc_map;
};




#endif //SOLVER_LINEARSOLVER_H
