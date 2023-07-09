//
// Created by gordan on 5/16/23.
//

#ifndef SOLVER_NEIGHBORPOSTPROCESSOR_H
#define SOLVER_NEIGHBORPOSTPROCESSOR_H

#include <variant>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/grid_in.h>
#include "CellPostprocessor.h"
#include "NuCurve.h"

using namespace dealii;

template <int dim>
class EggShellPostprocessor : public CellPostprocessor<dim> {

public:
    EggShellPostprocessor(unsigned int, unsigned int);
    void process(const Triangulation<dim>&  triangulation,
                 const Vector<double>&      solution,
                 const FE_Q<dim>&           fe,
                 std::vector<double>& result);

private:
    const Triangulation<dim> *triangulation_ptr = nullptr;
    const Vector<double> *solution_ptr = nullptr;
    const FE_Q<dim> *fe_ptr = nullptr;

    unsigned int egg_id;
    unsigned int shell_id;

};

#endif //SOLVER_NEIGHBORPOSTPROCESSOR_H