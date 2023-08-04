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

//
// Created by gordan on 6/15/23.
//

#ifndef SOLVER_CIRCLEPERIODICITYMAPPER_H
#define SOLVER_CIRCLEPERIODICITYMAPPER_H

#include <vector>
#include <limits>
#include <cmath>
#include <map>
#include <set>
#include <cassert>
#include <algorithm>
#include <cmath>
#include "IPeriodicityMapper.h"


using namespace std;

class CirclePeriodicityMapper : public IPeriodicityMapper{

public:
    CirclePeriodicityMapper(vector<unsigned int>,
                            vector<unsigned int>,
                            map<unsigned int, vector<double>>);

    set<pair<unsigned int, unsigned int>> get_matched_pair_indices() override;
    map<unsigned int, vector<double>> reference_to_center(map<unsigned int, vector<double>>);
    double get_weigth();
    void set_weigth(double);
    double error_overlap();

};


#endif //SOLVER_CIRCLEPERIODICITYMAPPER_H
