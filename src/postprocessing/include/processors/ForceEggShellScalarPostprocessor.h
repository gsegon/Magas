//
// Created by gordan on 5/16/23.
//

#ifndef SOLVER_FORCEEGGSHELLSCALARPOSTPROCESSOR_H
#define SOLVER_FORCEEGGSHELLSCALARPOSTPROCESSOR_H

#include <variant>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/grid_in.h>
#include "ScalarPostprocessor.h"
#include "NuCurve.h"

using namespace dealii;

template <int dim>
class ForceEggShellScalarPostprocessor : public ScalarPostprocessor<dim> {

public:
    ForceEggShellScalarPostprocessor(unsigned int, unsigned int, const std::string&);
    void process(const Triangulation<dim>&  triangulation,
                 const Vector<double>&      solution,
                 const FE_Q<dim>&           fe,
                 double& result);

private:
    const Triangulation<dim> *triangulation_ptr = nullptr;
    const Vector<double> *solution_ptr = nullptr;
    const FE_Q<dim> *fe_ptr = nullptr;

    unsigned int egg_id;
    unsigned int shell_id;
    std::string component;

};

#endif //SOLVER_FORCEEGGSHELLSCALARPOSTPROCESSOR_H