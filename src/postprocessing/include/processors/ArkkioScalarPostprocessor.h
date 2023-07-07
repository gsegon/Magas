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
#include "BHCurve.h"

using namespace dealii;

template <int dim>
class ArkkioScalarPostprocessor : public ScalarPostprocessor<dim> {

public:
    ArkkioScalarPostprocessor(unsigned int, const std::unordered_map<int, BHCurve*>&);

    void process(const Triangulation<dim>&  triangulation,
                 const Vector<double>&      solution,
                 const FE_Q<dim>&           fe,
                 double& result);

private:
    unsigned int mat_id{};
    const std::unordered_map<int, BHCurve*>* nu_map_ptr = nullptr;

};

#endif //SOLVER_ARKKIOSCALARPOSTPROCESSOR_H