//
// Created by gordan on 7/13/23.
//

#ifndef MAGAS_FSOURCE_H
#define MAGAS_FSOURCE_H

class FSource{

public:
    virtual double get_value(double x, double y) = 0;
};

#endif //MAGAS_FSOURCE_H
