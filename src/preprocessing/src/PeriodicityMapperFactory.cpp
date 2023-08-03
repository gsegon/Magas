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
// Created by gordan on 7/3/23.
//
#include <vector>
#include <regex>

#include "IPeriodicityMapper.h"
#include "PeriodicityMapperFactory.h"
#include "LinePeriodicityMapper.h"
#include "CirclePeriodicityMapper.h"


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