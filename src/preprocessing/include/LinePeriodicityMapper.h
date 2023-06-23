//
// Created by gordan on 6/15/23.
//

#ifndef SOLVER_LINEPERIODICITYMAPPER_H
#define SOLVER_LINEPERIODICITYMAPPER_H

#include <vector>
#include <limits>
#include <cmath>
#include <map>
#include "IPeriodicityMapper.h"


class LinePeriodicityMapper : public IPeriodicityMapper{

public:
    LinePeriodicityMapper(std::vector<unsigned int>,
                            std::vector<unsigned int>,
                            std::map<unsigned int, std::vector<double>>);

    std::set<std::pair<unsigned int, unsigned int>> get_matched_pair_indices() override;

};


LinePeriodicityMapper::LinePeriodicityMapper(   std::vector<unsigned int> a_dofs,
                                                std::vector<unsigned int> b_dofs,
                                                std::map<unsigned int, std::vector<double>> dof_to_nodes){
    this->a_dofs = a_dofs;
    this->b_dofs = b_dofs;
    this->dof_to_nodes = dof_to_nodes;

    assert(a_dofs.size() == b_dofs.size() && "Number of dofs 'a' not equal to number dofs 'b'.");

    std::vector<std::vector<double>> a_points;
    std::vector<std::vector<double>> b_points;

    std::vector<std::pair<unsigned int, double>> dofs_1;
    dofs_1.reserve(a_dofs.size());
    for(auto dof : a_dofs){
        dofs_1.emplace_back(dof, std::pow(dof_to_nodes[dof][0], 2) + std::pow(dof_to_nodes[dof][1], 2));
    }

    std::vector<std::pair<unsigned int, double>> dofs_2;
    dofs_2.reserve(b_dofs.size());
    for(auto dof : b_dofs){
        dofs_2.emplace_back(dof, std::pow(dof_to_nodes[dof][0], 2) + std::pow(dof_to_nodes[dof][1], 2));
    }

    std::sort(dofs_1.begin(), dofs_1.end(), [](std::pair<unsigned int, double> a, std::pair<unsigned int, double> b) {return std::get<1>(a) < std::get<1>(b);});
    std::sort(dofs_2.begin(), dofs_2.end(), [](std::pair<unsigned int, double> a, std::pair<unsigned int, double> b) {return std::get<1>(a) < std::get<1>(b);});

    for (int i =0; i < (int)dofs_1.size(); i++) {
        auto pair = std::pair<unsigned int, unsigned int>(std::get<0>(dofs_1[i]), std::get<0>(dofs_2[i]));
        matched_pairs.insert(pair);
    }

}

std::set<std::pair<unsigned int, unsigned int>> LinePeriodicityMapper::get_matched_pair_indices() {
    return this->matched_pairs;
}


#endif //SOLVER_LINEPERIODICITYMAPPER_H
