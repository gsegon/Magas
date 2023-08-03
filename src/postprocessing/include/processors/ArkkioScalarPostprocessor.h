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


#ifndef SOLVER_ARKKIOSCALARPOSTPROCESSOR_H
#define SOLVER_ARKKIOSCALARPOSTPROCESSOR_H

#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/grid_in.h>
#include "ScalarPostprocessor.h"
#include "NuCurve.h"

using namespace dealii;

template <int dim>
class ArkkioScalarPostprocessor : public ScalarPostprocessor<dim> {

public:
    ArkkioScalarPostprocessor(unsigned int, const std::unordered_map<int, NuCurve*>&);

    void process(const Triangulation<dim>&  triangulation,
                 const Vector<double>&      solution,
                 const FE_Q<dim>&           fe,
                 double& result);

private:
    unsigned int mat_id{};
    const std::unordered_map<int, NuCurve*>* nu_map_ptr = nullptr;

};

#endif //SOLVER_ARKKIOSCALARPOSTPROCESSOR_H