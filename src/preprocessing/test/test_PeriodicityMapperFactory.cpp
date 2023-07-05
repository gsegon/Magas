//
// Created by gordan on 3/12/23.
//
#include <gtest/gtest.h>

#include "PeriodicityMapperFactory.h"

TEST(TestPeriodicityMapperFactory, basic) {

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

    PeriodicityMapperFactory pmf{a_dofs, b_dofs, dof_to_nodes};
    auto a = pmf.create("periodic-line");
    auto pairs = a->get_matched_pair_indices();

    ASSERT_EQ(pairs, truth);

}

