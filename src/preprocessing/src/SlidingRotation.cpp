//
// Created by gordan on 7/11/23.
//

#include <algorithm>
#include <cmath>

#include "SlidingRotation.h"

SlidingRotation::SlidingRotation(vector<dof> dofs, map<dof, vector<double>> dof_to_node, int offset) {
    this->dofs = dofs;
    this->dof_to_node = dof_to_node;
    this->offset = offset;

    this->sort();
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

vector<dof> SlidingRotation::get_dofs() {
    return dofs;
}

dof SlidingRotation::get_mapped(dof given) {
    auto pos_given = std::find(dofs.begin(), dofs.end(), given) - dofs.begin();

    auto pos_mapped = pos_given + offset;
    pos_mapped = pos_mapped % dofs.size();

    if (pos_given >= dofs.size())
        return given;
    else
        return dofs[pos_mapped];

}