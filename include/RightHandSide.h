//
// Created by gordan on 5/11/23.
//

#ifndef SOLVER_RIGHTHANDSIDE_H
#define SOLVER_RIGHTHANDSIDE_H


class RightHandSide{
public:
    virtual double value(int physical_id);
};



#endif //SOLVER_RIGHTHANDSIDE_H
