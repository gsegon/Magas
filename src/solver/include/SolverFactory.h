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

//
// Created by gordan on 7/14/23.
//

#ifndef MAGAS_SOLVERFACTORY_H
#define MAGAS_SOLVERFACTORY_H

#include "Solver.h"
#include "NewtonSolver.h"
#include "LinearSolver.h"


template<int dim>
class SolverFactory {

public:
   Solver<dim>* create(bool nonlinear){
       if (nonlinear) return new NewtonSolver<dim>();
       else return new LinearSolver<dim>();
    }
};


#endif //MAGAS_SOLVERFACTORY_H
