//
// Created by gordan on 7/13/23.
//

#ifndef MAGAS_CONSTFSOURCE_H
#define MAGAS_CONSTFSOURCE_H

#include "FSource.h"

class ConstFSource : public FSource {

public:
    ConstFSource(double value);
    double get_value(double x, double y);

private:
    double value;

};


#endif //MAGAS_CONSTFSOURCE_H
