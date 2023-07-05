//
// Created by gordan on 7/3/23.
//
#include <regex>

#include "IPeriodicityMapper.h"
#include "PeriodicityMapperFactory.h"
#include "LinePeriodicityMapper.h"
#include "CirclePeriodicityMapper.h"
#include <vector>

using namespace std;

PeriodicityMapperFactory::PeriodicityMapperFactory(vector<unsigned int> a_dofs,
                                                   vector<unsigned int> b_dofs,
                                                   map<unsigned int, vector<double>> dof_to_nodes){
    this->a_dofs = a_dofs;
    this->b_dofs = b_dofs;
    this->dof_to_nodes = dof_to_nodes;
}


IPeriodicityMapper* PeriodicityMapperFactory::create(std::string input_string) {

    if (input_string == "periodic-line"){
        auto mapper = new LinePeriodicityMapper{a_dofs, b_dofs, dof_to_nodes};
        mapper->set_weigth(1);
        return mapper;
    }
    if (input_string == "anti-periodic-line"){
        auto mapper = new LinePeriodicityMapper{a_dofs, b_dofs, dof_to_nodes};
        mapper->set_weigth(-1);
        return mapper;
    }
    if (input_string == "periodic-circle"){
        auto mapper = new CirclePeriodicityMapper{a_dofs, b_dofs, dof_to_nodes};
        mapper->set_weigth(1);
        return mapper;
    }
    if (input_string == "anti-periodic-circle"){
        auto mapper = new CirclePeriodicityMapper{a_dofs, b_dofs, dof_to_nodes};
        mapper->set_weigth(-1);
        return mapper;
    }

    runtime_error("PeriodicityMapperFactory cannot create " + input_string + " type mapper.");
    return nullptr;
}