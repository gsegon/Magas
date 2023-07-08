//
// Created by gordan on 7/3/23.
//

#ifndef MAGAS_NUCURVEFACTORY_H
#define MAGAS_NUCURVEFACTORY_H

#include <map>
#include <variant>
#include <vector>
#include <filesystem>

#include "NuCurve.h"

using namespace std;

class NuCurveFactory {

public:
    NuCurve* create(double input);
    NuCurve* create(string input);
    NuCurve* create(filesystem::path input);

};


#endif //MAGAS_NUCURVEFACTORY_H
