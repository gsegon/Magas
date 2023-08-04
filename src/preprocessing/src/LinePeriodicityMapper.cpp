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
// Created by gordan on 7/5/23.
//

#include <map>
#include <variant>
#include <vector>
#include <algorithm>
#include <stdexcept>

#include "LinePeriodicityMapper.h"

LinePeriodicityMapper::LinePeriodicityMapper(   vector<unsigned int> a_dofs,
                                                vector<unsigned int> b_dofs,
                                                map<unsigned int, vector<double>> dof_to_nodes){
    this->a_dofs = a_dofs;
    this->b_dofs = b_dofs;
    this->dof_to_nodes = dof_to_nodes;

    assert(a_dofs.size() == b_dofs.size() && "Number of dofs 'a' not equal to number dofs 'b'.");

    vector<vector<double>> a_points;
    vector<vector<double>> b_points;

    vector<pair<unsigned int, double>> dofs_1;
    dofs_1.reserve(a_dofs.size());
    for(auto dof : a_dofs){
        dofs_1.emplace_back(dof, pow(dof_to_nodes[dof][0], 2) + pow(dof_to_nodes[dof][1], 2));
    }

    vector<pair<unsigned int, double>> dofs_2;
    dofs_2.reserve(b_dofs.size());
    for(auto dof : b_dofs){
        dofs_2.emplace_back(dof, pow(dof_to_nodes[dof][0], 2) + pow(dof_to_nodes[dof][1], 2));
    }

    sort(dofs_1.begin(), dofs_1.end(), [](pair<unsigned int, double> a, pair<unsigned int, double> b) {return get<1>(a) < get<1>(b);});
    sort(dofs_2.begin(), dofs_2.end(), [](pair<unsigned int, double> a, pair<unsigned int, double> b) {return get<1>(a) < get<1>(b);});

    for (int i =0; i < (int)dofs_1.size(); i++) {
        auto math_pair = pair<unsigned int, unsigned int>(get<0>(dofs_1[i]), get<0>(dofs_2[i]));
        matched_pairs.insert(math_pair);
    }

    for (auto pair: matched_pairs){
        auto r1 = sqrt(pow(dof_to_nodes[pair.first][0], 2) + pow(dof_to_nodes[pair.first][1], 2));
        auto r2 = sqrt(pow(dof_to_nodes[pair.second][0], 2) + pow(dof_to_nodes[pair.second][1], 2));
        auto e_distance = abs(r1-r2);
        if (e_distance > 1e-8){
            throw std::runtime_error("Radius of matched points too big.");
        }
    }

}

set<pair<unsigned int, unsigned int>> LinePeriodicityMapper::get_matched_pair_indices() {
    return this->matched_pairs;
}

double LinePeriodicityMapper::get_weigth() {
    return weigth;
}

void LinePeriodicityMapper::set_weigth(double val){
    weigth = val;
}