//
// Created by gordan on 5/16/23.
//

#ifndef SOLVER_FLUXLINKAGESCALARPOSTPROCESSOR_H
#define SOLVER_FLUXLINKAGESCALARPOSTPROCESSOR_H

#include <variant>

#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/grid_in.h>

#include "processors/ScalarPostprocessor.h"
#include "FSource.h"

using namespace dealii;

template <int dim>
class FluxLinkageScalarPostprocessor : public ScalarPostprocessor<dim> {

public:
    FluxLinkageScalarPostprocessor(unsigned int, std::unordered_map<int, std::variant<FSource*, std::pair<double, double>>>&);

    void process(const Triangulation<dim>&  triangulation,
                 const Vector<double>&      solution,
                 const FE_Q<dim>&           fe,
                 double& result);

private:
    unsigned int mat_id{};
    std::unordered_map<int, std::variant<FSource*, std::pair<double, double>>>* f_map_ptr = nullptr;

};

#endif //SOLVER_FLUXLINKAGESCALARPOSTPROCESSOR_H