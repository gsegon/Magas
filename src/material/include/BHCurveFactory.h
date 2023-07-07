//
// Created by gordan on 7/3/23.
//

#ifndef MAGAS_BHCURVEFACTORY_H
#define MAGAS_BHCURVEFACTORY_H

#include <map>
#include <variant>
#include <vector>

#include "BHCurve.h"

using namespace std;

class BHCurveFactory {

public:
    BHCurve* create(double input);
    BHCurve* create(string input);

};


#endif //MAGAS_BHCURVEFACTORY_H
