//
// Created by gordan on 7/3/23.
//

#ifndef MAGAS_PERIODICIYMAPPERFACTORY_H
#define MAGAS_PERIODICIYMAPPERFACTORY_H

#include <map>
#include <variant>
#include <vector>

#include "IPeriodicityMapper.h"

//using namespace std;

class PeriodicityMapperFactory {

public:
    PeriodicityMapperFactory(std::vector<unsigned int> a_dofs,
                             std::vector<unsigned int> b_dofs,
                             std::map<unsigned int, std::vector<double>> dof_to_nodes);
    IPeriodicityMapper* create(std::string input_string);

private:
    std::vector<unsigned int> a_dofs;
    std::vector<unsigned int> b_dofs;
    std::map<unsigned int, std::vector<double>> dof_to_nodes;

};


#endif //MAGAS_PERIODICIYMAPPERFACTORY_H
