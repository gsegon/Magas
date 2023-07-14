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
