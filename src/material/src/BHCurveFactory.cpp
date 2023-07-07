//
// Created by gordan on 7/3/23.
//
#include <regex>

#include <map>
#include <variant>
#include <vector>

#include "../include/BHCurveFactory.h"
#include "../include/LinearBHCurve.h"
#include "../include/AnalyticBHCurve.h"
#include "../include/BHCurve.h"

using namespace std;


BHCurve* BHCurveFactory::create(double input) {
    return new LinearBHCurve{input};
}

BHCurve* BHCurveFactory::create(string input) {
    if (input == "Analytic") return new AnalyticBHCurve{};
    throw runtime_error("No BHCurve can be created for \'" + input + "\' input.");
}