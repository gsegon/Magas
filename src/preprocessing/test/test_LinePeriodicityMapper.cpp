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
// Created by gordan on 3/12/23.
//
#include <gtest/gtest.h>

#include "LinePeriodicityMapper.h"

TEST(TestLinePeriodicityMapper, basic) {

    std::set<std::pair<unsigned int, unsigned int>> truth{
                                                              {1, 5},
                                                              {0, 4},
                                                              {2, 7},
                                                              {3, 6}};

    std::map<unsigned int, std::vector<double>> dof_to_nodes;

    dof_to_nodes[0] = {1.0, 0.0};
    dof_to_nodes[1] = {1.5, 0.0};
    dof_to_nodes[2] = {2.0, 0.0};
    dof_to_nodes[3] = {1.8, 0.0};

    dof_to_nodes[4] = {0.0, 1.0};
    dof_to_nodes[5] = {0.0, 1.5};
    dof_to_nodes[6] = {0.0, 1.8};
    dof_to_nodes[7] = {0.0, 2.0};

    std::vector<unsigned int> a_dofs{1, 0, 2, 3};
    std::vector<unsigned int> b_dofs{4, 5, 6, 7};

    LinePeriodicityMapper lpm{a_dofs, b_dofs, dof_to_nodes};
    auto pairs = lpm.get_matched_pair_indices();

    ASSERT_EQ(pairs, truth);

}

TEST(TestLinePeriodicityMapper, distance_too_big) {

    std::set<std::pair<unsigned int, unsigned int>> truth{
            {1, 5},
            {0, 4},
            {2, 7},
            {3, 6}};

    std::map<unsigned int, std::vector<double>> dof_to_nodes;

    dof_to_nodes[0] = {1.0, 0.0};
    dof_to_nodes[1] = {1.5, 0.0};
    dof_to_nodes[2] = {2.0, 0.0};
    dof_to_nodes[3] = {1.8, 0.0};

    double eps = 1e-3;
    dof_to_nodes[4] = {0.0, 1.0};
    dof_to_nodes[5] = {0.0, 1.5};
    dof_to_nodes[6] = {0.0, 1.8+eps};
    dof_to_nodes[7] = {0.0, 2.0};

    std::vector<unsigned int> a_dofs{1, 0, 2, 3};
    std::vector<unsigned int> b_dofs{4, 5, 6, 7};

    try{
        LinePeriodicityMapper lpm{a_dofs, b_dofs, dof_to_nodes};
        auto pairs = lpm.get_matched_pair_indices();
        ASSERT_TRUE(0);
    }
    catch(const runtime_error& error){
        cout << "Exception thrown (expected): " << error.what() << endl;
        ASSERT_TRUE(1);
    }

}
