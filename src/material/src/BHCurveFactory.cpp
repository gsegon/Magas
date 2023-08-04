// Magas - Magnetostatic Analysis Suite
// Copyright (C) 2023  Gordan Segon
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
// USA

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