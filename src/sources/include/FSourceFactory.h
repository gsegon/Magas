//
// Created by gordan on 7/13/23.
//

#ifndef MAGAS_FSOURCEFACTORY_H
#define MAGAS_FSOURCEFACTORY_H

#include <string>
#include "FSource.h"

class FSourceFactory {

public:
    FSource* create(double val);
    FSource* create(std::string);

};


#endif //MAGAS_FSOURCEFACTORY_H
