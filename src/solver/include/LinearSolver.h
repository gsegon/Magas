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
#include "BHCurve.h"

using namespace dealii;

template <int dim>
class LinearSolver{

public:
    LinearSolver();
    void read_mesh(std::string);
    void setup_system();
    void assemble_system();
    void solve();
    void set_nu_map(std::unordered_map<int, BHCurve*>);
    void set_f_map(std::unordered_map<int, std::variant<double, std::pair<double, double>>>);
    void set_dc_map(std::unordered_map<int, double>);
    void set_per_map(std::unordered_map<std::string, std::vector<unsigned int>>);

    Triangulation<dim>& get_triangulation();
    Vector<double>& get_solution();
    Vector<double>& get_rhs();
    FE_Q<dim>& get_fe();

private:

    Triangulation<dim> triangulation;
    FE_Q<dim> fe;
    DoFHandler<dim> dof_handler;
    QGauss<dim> quadrature_formula;

    AffineConstraints<double> constraints;

    struct AssemblyScratchData{
        AssemblyScratchData(const FiniteElement<dim> &fe);
        AssemblyScratchData(const AssemblyScratchData &scratch_data);

        FEValues<dim> fe_values;
        std::vector<double> rhs_values;
    };

    struct AssemblyCopyData{
        FullMatrix<double> cell_matrix;
        Vector<double> cell_rhs;
        std::vector<types::global_dof_index> local_dof_indices;
    };

    void local_assemble_system(const typename DoFHandler<dim>::active_cell_iterator& cell,
                               AssemblyScratchData& scratch,
                               AssemblyCopyData& copy_data);
    void copy_local_to_global(const AssemblyCopyData &copy_data);

    SparsityPattern sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double> solution;
    Vector<double> system_rhs;
    std::unordered_map<int, BHCurve*> nu_map;
    std::unordered_map<int, std::variant<double, std::pair<double, double>>> f_map;
    std::unordered_map<int, double> dc_map;
    std::unordered_map<std::string, std::vector<unsigned int>> per_map;
};

#endif //SOLVER_LINEARSOLVER_H
