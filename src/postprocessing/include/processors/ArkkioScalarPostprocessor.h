//
// Created by gordan on 5/16/23.
//

#ifndef SOLVER_ARKKIOSCALARPOSTPROCESSOR_H
#define SOLVER_ARKKIOSCALARPOSTPROCESSOR_H

#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/grid_in.h>
#include "ScalarPostprocessor.h"

using namespace dealii;

template <int dim>
class ArkkioScalarPostprocessor : public ScalarPostprocessor<dim> {

public:
    ArkkioScalarPostprocessor(unsigned int, const std::unordered_map<int, double>&);

    void process(const Triangulation<dim>&  triangulation,
                 const Vector<double>&      solution,
                 const FE_Q<dim>&           fe,
                 double& result);

private:
//    const Triangulation<dim> *triangulation_ptr = nullptr;
//    const Vector<double> *solution_ptr = nullptr;
//    const FE_Q<dim> *fe_ptr = nullptr;

    unsigned int mat_id{};
    const std::unordered_map<int, double>* nu_map_ptr = nullptr;

};

#endif //SOLVER_ARKKIOSCALARPOSTPROCESSOR_H