//
// Created by gordan on 7/3/23.
//
#include <regex>

#include <map>
#include <variant>
#include <vector>
#include <filesystem>

#include "NuCurveFactory.h"
#include "LinearNuCurve.h"
#include "AnalyticNuCurve.h"
#include "InterpolatedNuCurve.h"
#include "NuCurve.h"
#include "strtk.hpp"

using namespace std;


NuCurve* NuCurveFactory::create(double input) {
    return new LinearNuCurve{input};
}

NuCurve* NuCurveFactory::create(string input) {
    if (input == "Analytic") return new AnalyticNuCurve{};
    throw runtime_error("No BHCurve can be created for \'" + input + "\' input.");
}

NuCurve* NuCurveFactory::create(filesystem::path input) {

    std::vector<double> b;
    std::vector<double> h;

    strtk::token_grid::options options;
    options.column_delimiters =",";

    strtk::token_grid bh_grid(input, options);
    for (std::size_t r = 0; r < bh_grid.row_count(); ++r){
        strtk::token_grid::row_type row = bh_grid.row(r);

        if (r==0){
        }
        else{
            b.push_back(std::stod(row.get<std::string>(0)));
            h.push_back(std::stod(row.get<std::string>(1)));
        }
        std::cout << std::endl;
    }
    return new InterpolatedNuCurve(b, h);
}