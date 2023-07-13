//
// Created by gordan on 7/13/23.
//

#ifndef MAGAS_HSOURCE_H
#define MAGAS_HSOURCE_H

#include <vector>

class HSource{

    virtual std::vector<double> get_value(double x, double y) = 0;

};

#endif //MAGAS_HSOURCE_H
