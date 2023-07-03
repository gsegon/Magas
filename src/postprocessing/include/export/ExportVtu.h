//
// Created by gordan on 5/15/23.
//

#ifndef SOLVER_EXPORTVTU_H
#define SOLVER_EXPORTVTU_H

#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/grid_in.h>
#include "processors/CellPostprocessor.h"

using namespace dealii;

template <int dim>
class ExportVtu {

public:
    ExportVtu(const Triangulation<dim>& triangulation);
    ExportVtu(const Triangulation<dim>& triangulation, const Vector<double>& rhs, const Vector<double>& solution, const FE_Q<dim>& fe);

    void attach_postprocessor(CellPostprocessor<dim>*, std::string);

    void write(const std::string&);


protected:
    const Triangulation<dim>* triangulation_ptr = nullptr;
    const Vector<double>* solution_ptr = nullptr;
    const Vector<double>* rhs_ptr = nullptr;
    const FE_Q<dim>* fe_ptr = nullptr;
    std::unordered_map<std::string, CellPostprocessor<dim>*> post_processors;
};


#endif //SOLVER_EXPORTVTU_H
