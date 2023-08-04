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
