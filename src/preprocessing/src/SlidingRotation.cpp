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
// Created by gordan on 7/11/23.
//

#include <iostream>
#include <algorithm>
#include <cmath>

#include "SlidingRotation.h"

SlidingRotation::SlidingRotation(vector<dof> dofs, map<dof, vector<double>> dof_to_node, int offset) {
    this->dofs = dofs;
    this->dof_to_node = dof_to_node;
    this->offset = offset;

    this->sort();
}

bool SlidingRotation::is_full_circle(){
    double x=0;
    double y=0;
    for (auto [dof, node] : dof_to_node){
        x += node[0];
        y += node[1];
    }
    x /= dofs.size();
    y /= dofs.size();

    if (std::sqrt(std::pow(x,2) + std::pow(y, 2)) < 1e-6)
        return true;
    else
        return false;

}

bool SlidingRotation::is_in(dof given) {
    return std::find(dofs.begin(), dofs.end(), given) - dofs.begin() < dofs.size();
}

void SlidingRotation::sort() {

    std::map<unsigned int, std::vector<double>> circle_points;
    for(auto dof : dofs)
        circle_points[dof] = dof_to_node[dof];

    std::vector<std::pair<unsigned int, double>> thetas;
    thetas.reserve(circle_points.size());
    for(auto [dof, point] : circle_points){
        thetas.emplace_back(dof, std::atan2(point[1], point[0]));
    }

    std::sort(thetas.begin(), thetas.end(), [](std::pair<unsigned int, double> a, std::pair<unsigned int, double> b) {return std::get<1>(a) < std::get<1>(b);});

    vector<dof> new_dofs;
    for (auto [dof, theta] : thetas){
        new_dofs.push_back(dof);
    }
    dofs = new_dofs;
}

void SlidingRotation::print_map() {
    std::cout << "Printout of Sliding rotation map: "<< std::endl;
    for (auto dof : dofs){
        std::cout << "Dof " << dof << "->" << get_mapped(dof) << std::endl;
    }
}

vector<dof> SlidingRotation::get_dofs() {
    return dofs;
}


dof SlidingRotation::get_mapped(dof given) {
    auto pos_given = std::find(dofs.begin(), dofs.end(), given) - dofs.begin();
    if (pos_given >= dofs.size())
        return given;

    auto pos_mapped = pos_given + offset;

    if (!is_full_circle()) {
        pos_mapped = pos_mapped % (dofs.size()-1);
    }
    else
        pos_mapped = pos_mapped % (dofs.size());

    return dofs[pos_mapped];

}