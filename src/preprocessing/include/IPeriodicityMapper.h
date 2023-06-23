//
// Created by gordan on 6/23/23.
//

#ifndef SOLVER_IPERIODICITYMAPPER_H
#define SOLVER_IPERIODICITYMAPPER_H


#include <vector>

class IPeriodicityMapper{

public:
//    IPeriodicityMapper(std::vector<unsigned int>, std::vector<unsigned int>, std::vector<std::vector<double>>);
    virtual std::set<std::pair<unsigned int, unsigned int>> get_matched_pair_indices() = 0;

protected:
    std::vector<unsigned int> a_dofs;
    std::vector<unsigned int> b_dofs;
    std::map<unsigned int, std::vector<double>> dof_to_nodes;
    std::set<std::pair<unsigned int, unsigned int>> matched_pairs;

};


#endif //SOLVER_IPERIODICITYMAPPER_H
