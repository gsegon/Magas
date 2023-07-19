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
#include "SlidingRotation.h"

using namespace dealii;

template <int dim>
class LinearSolver: public Solver<dim>{

public:
    LinearSolver();
    void read_mesh(const std::string&);
    void setup_system();
    void setup_rotation(unsigned int, unsigned int, int);
    void extend_dsp(DynamicSparsityPattern& dsp);
    void assemble_system();
    void solve();
    void set_nu_map(t_nu_map);
    void set_f_map(t_f_map);
    void set_dc_map(t_dc_map);
    void set_per_map(t_per_map);
    void set_rot_map(t_rot_map);
    void run();

    Triangulation<dim>& get_triangulation();
    Vector<double>& get_solution();
    Vector<double>& get_rhs();
    FE_Q<dim>& get_fe();

private:

    std::vector<unsigned int> rot_cell_indices;
    std::vector<unsigned int> rot_dofs;
    SlidingRotation* sr;

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
    t_nu_map nu_map;
    t_f_map f_map;
    t_dc_map dc_map;
    t_per_map per_map;
    t_rot_map rot_map;

};

#endif //SOLVER_LINEARSOLVER_H