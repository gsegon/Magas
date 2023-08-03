// Magas - Magnetostatic Analysis Suite
// Copyright (C) 2023  Gordan Segon
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
// USA

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
