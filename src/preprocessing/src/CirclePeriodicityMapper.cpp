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

#include "CirclePeriodicityMapper.h"
#include <stdexcept>


map<unsigned int, vector<double>> CirclePeriodicityMapper::reference_to_center(map<unsigned int, vector<double>> points){
    vector<double> centeroid{0,0};
    for (auto [dof, point] : points){
        centeroid[0] += point[0];
        centeroid[1] += point[1];
    }
    centeroid[0] /= points.size();
    centeroid[1] /= points.size();

    map<unsigned int, vector<double>> new_points;
    for (auto [dof, point] : points){
        new_points[dof] = {point[0]-centeroid[0], point[1]-centeroid[1]};
    }

    return new_points;
}


CirclePeriodicityMapper::CirclePeriodicityMapper(vector<unsigned int> a_dofs,
                                                 vector<unsigned int> b_dofs,
                                                 map<unsigned int, vector<double>> dof_to_nodes){
    this->a_dofs = a_dofs;
    this->b_dofs = b_dofs;
    this->dof_to_nodes = dof_to_nodes;

    assert(a_dofs.size() == b_dofs.size() && "Number of dofs 'a' not equal to number dofs 'b'.");

    map<unsigned int, vector<double>> a_points;
    map<unsigned int, vector<double>> b_points;

    for(auto dof : a_dofs)
        a_points[dof] = dof_to_nodes[dof];

    for(auto dof : b_dofs)
        b_points[dof] = dof_to_nodes[dof];

    a_points = reference_to_center(a_points);
    b_points = reference_to_center(b_points);

    std::vector<std::pair<unsigned int, double>> thetas_1;
    thetas_1.reserve(a_dofs.size());
    for(auto [dof, point] : a_points){
        thetas_1.emplace_back(dof, atan2(point[1], point[0]));
    }

    std::vector<std::pair<unsigned int, double>> thetas_2;
    thetas_2.reserve(a_dofs.size());
    for(auto [dof, point] : b_points){
        thetas_2.emplace_back(dof, atan2(point[1], point[0]));
    }

    std::sort(thetas_1.begin(), thetas_1.end(), [](std::pair<unsigned int, double> a, std::pair<unsigned int, double> b) {return std::get<1>(a) < std::get<1>(b);});
    std::sort(thetas_2.begin(), thetas_2.end(), [](std::pair<unsigned int, double> a, std::pair<unsigned int, double> b) {return std::get<1>(a) < std::get<1>(b);});

    for (int i =0; i < (int)thetas_1.size(); i++) {
        auto pair = std::pair<unsigned int, unsigned int>(std::get<0>(thetas_1[i]), std::get<0>(thetas_2[i]));
        matched_pairs.insert(pair);
    }

    for (auto pair: matched_pairs){
        auto e_distance = sqrt(pow((a_points[pair.first][0]-b_points[pair.second][0]), 2) + pow((a_points[pair.first][1]-b_points[pair.second][1]), 2));
        if (e_distance > 1e-8){
            throw std::runtime_error("Distance of matched points too big.");
        }
    }

}


set<pair<unsigned int, unsigned int>> CirclePeriodicityMapper::get_matched_pair_indices() {
    return this->matched_pairs;
}

double CirclePeriodicityMapper::get_weigth() {
    return weigth;
}

void CirclePeriodicityMapper::set_weigth(double val){
    weigth = val;
}