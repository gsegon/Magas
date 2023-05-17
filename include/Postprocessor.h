//
// Created by gordan on 5/16/23.
//

#ifndef SOLVER_POSTPROCESSOR_H
#define SOLVER_POSTPROCESSOR_H

#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/grid_in.h>

using namespace dealii;

template <int dim>
class Postprocessor {

    public:
        virtual void process(const Triangulation<dim> &triangulation,
                             const Vector<double>& solution,
                             const FE_Q<dim>& fe,
                             std::vector<double>& result) = 0;

};

#endif //SOLVER_POSTPROCESSOR_H