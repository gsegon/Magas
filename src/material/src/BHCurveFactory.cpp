//
// Created by gordan on 7/3/23.
//
#include <regex>

#include <map>
#include <variant>
#include <vector>

#include "../include/NuCurveFactory.h"
#include "../include/LinearNuCurve.h"
#include "../include/AnalyticNuCurve.h"
#include "../include/NuCurve.h"

using namespace std;


NuCurve* NuCurveFactory::create(double input) {
    return new LinearNuCurve{input};
}

NuCurve* NuCurveFactory::create(string input) {
    if (input == "Analytic") return new AnalyticNuCurve{};
    throw runtime_error("No BHCurve can be created for \'" + input + "\' input.");
}