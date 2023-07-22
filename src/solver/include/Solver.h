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
#include <deal.II/numerics/data_out.h>

#include "NuCurveFactory.h"
#include "FSourceFactory.h"
#include "SlidingRotation.h"

using namespace dealii;

typedef std::unordered_map<int, NuCurve*> t_nu_map;
typedef std::unordered_map<int, std::variant<FSource*, std::pair<double, double>>> t_f_map;
typedef std::unordered_map<int, double> t_dc_map;
typedef std::unordered_map<std::string, std::vector<unsigned int>> t_per_map;
typedef std::map<std::string, std::string> t_postprocessor_strings;

template<int dim>
class Solver{
public:

    Solver();
    void read_mesh(const std::string&);
    void set_nu_map(t_nu_map);
    void set_f_map(t_f_map);
    void set_dc_map(t_dc_map);
    void set_per_map(t_per_map);
    virtual void setup_system() = 0;
    void assemble_system();
    virtual void run() = 0;
    virtual void solve() = 0;

    Triangulation<dim>& get_triangulation();
    Vector<double>& get_solution();
    Vector<double>& get_rhs();
    FE_Q<dim>& get_fe();

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

    virtual void local_assemble_system(const typename DoFHandler<dim>::active_cell_iterator& cell,
                               AssemblyScratchData& scratch,
                               AssemblyCopyData& copy_data) = 0;
    void copy_local_to_global(const AssemblyCopyData &copy_data);

    SparsityPattern sparsity_pattern;
    SparseMatrix<double> system_matrix;

    std::vector<unsigned int> rot_cell_indices;
    std::vector<unsigned int> rot_dofs;
    SlidingRotation* sr = nullptr;

    Vector<double> solution;
    Vector<double> system_rhs;
    t_nu_map nu_map;
    t_f_map f_map;
    t_dc_map dc_map;
    t_per_map per_map;
    t_rot_map rot_map = {};
};

#endif //MAGAS_SOLVER_H