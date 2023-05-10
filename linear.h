//
// Created by gordan on 5/10/23.
//

#ifndef SOLVER_LINEAR_H
#define SOLVER_LINEAR_H

#include <deal.II/grid/tria.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>


using namespace dealii;

template <int dim>
class LinearSolver{

public:
    LinearSolver();
    void run();

private:
    void read_mesh(std::string);
    void setup_system();
    void assemble_system();
    void solve();
    void output_results() const;

    Triangulation<dim> triangulation;
    FE_Q<dim> fe;
    DoFHandler<dim> dof_handler;

    SparsityPattern sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double> solution;
    Vector<double> system_rhs;
};

class RightHandSide{
public:
    virtual double value(int physical_id);
};



#endif //SOLVER_LINEAR_H
