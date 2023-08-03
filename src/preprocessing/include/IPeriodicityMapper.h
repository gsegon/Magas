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
// Created by gordan on 6/23/23.
//

#ifndef SOLVER_IPERIODICITYMAPPER_H
#define SOLVER_IPERIODICITYMAPPER_H


#include <vector>
#include <set>

class IPeriodicityMapper{

public:
    virtual std::set<std::pair<unsigned int, unsigned int>> get_matched_pair_indices() = 0;
    virtual double get_weigth() = 0;

protected:
    std::vector<unsigned int> a_dofs;
    std::vector<unsigned int> b_dofs;
    std::map<unsigned int, std::vector<double>> dof_to_nodes;
    std::set<std::pair<unsigned int, unsigned int>> matched_pairs;
    double weigth;

};


#endif //SOLVER_IPERIODICITYMAPPER_H
