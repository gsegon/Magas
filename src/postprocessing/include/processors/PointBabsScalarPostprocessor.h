//
// Created by gordan on 5/16/23.
//

#ifndef SOLVER_POINTBABSSCALARPOSTPROCESSOR_H
#define SOLVER_POINTBABSSCALARPOSTPROCESSOR_H

#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/grid_in.h>

using namespace dealii;

template <int dim>
class PointBabsScalarPostprocessor {

public:
    PointBabsScalarPostprocessor(std::vector<double>);

    void process(const Triangulation<dim>&  triangulation,
                 const Vector<double>&      solution,
                 const FE_Q<dim>&           fe,
                 double& result);

private:
    const Triangulation<dim> *triangulation_ptr = nullptr;
    const Vector<double> *solution_ptr = nullptr;
    const FE_Q<dim> *fe_ptr = nullptr;

    std::vector<double> point;

};

#endif //SOLVER_POINTBABSSCALARPOSTPROCESSOR_H