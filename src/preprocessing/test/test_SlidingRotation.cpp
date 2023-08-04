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

#include <gtest/gtest.h>
#include <tuple>
#include <unordered_map>
#include <string>
#include <any>
#include <fstream>
#include <filesystem>
#include <cmath>

#include "SlidingRotation.h"


TEST(TestSlidingRotation, basic) {

    std::map<unsigned int, std::vector<double>> dof_to_nodes;

    dof_to_nodes[1] = {1.0, 0.0};
    dof_to_nodes[2] = {0.0, 1.0};
    dof_to_nodes[3] = {-1.0, 0.0};
    dof_to_nodes[4] = {0, -1.0};

    std::vector<unsigned int> dofs = {1, 2, 3, 4};

    SlidingRotation sr{dofs, dof_to_nodes, 2};
    for (auto dof : sr.get_dofs()){
        std::cout << "dof: " << dof << " mapped: " << sr.get_mapped(dof) << std::endl;
    }
}

TEST(TestSlidingRotation, periodic_section) {

    std::map<unsigned int, std::vector<double>> dof_to_nodes;
    std::vector<unsigned int> dofs;

    for (int i=0; i<50; i++){
        dof_to_nodes[i]={cos(M_PI/2.0 * i/50.0), sin(M_PI/2.0 * i/50.0)};
        dofs.push_back(i);
    }

    SlidingRotation sr{dofs, dof_to_nodes, 2};
    sr.print_map();

}