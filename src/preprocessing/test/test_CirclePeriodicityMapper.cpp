//
// Created by gordan on 3/12/23.
//
#include <gtest/gtest.h>

#include "CirclePeriodicityMapper.h"

TEST(CirclePeriodicityMapper, translated) {

    std::set<std::pair<unsigned int, unsigned int>> truth{
                                                              {1, 5},
                                                              {0, 4},
                                                              {2, 6},
                                                              {3, 7}};

    std::map<unsigned int, std::vector<double>> dof_to_nodes;

    dof_to_nodes[0] = {1.0, 0.0};
    dof_to_nodes[1] = {0.0, 1.0};
    dof_to_nodes[2] = {-1.0, 0.0};
    dof_to_nodes[3] = {0, -1.0};

    double x_offset = 3.5;
    double y_offset = -0.5;
    dof_to_nodes[4] = {1.0+x_offset, 0.0+y_offset};
    dof_to_nodes[5] = {0.0+x_offset, 1.0+y_offset};
    dof_to_nodes[6] = {-1.0+x_offset, 0.0+y_offset};
    dof_to_nodes[7] = {0.0+x_offset, -1.0+y_offset};

    std::vector<unsigned int> a_dofs{1, 0, 2, 3};
    std::vector<unsigned int> b_dofs{4, 5, 6, 7};

    CirclePeriodicityMapper lpm{a_dofs, b_dofs, dof_to_nodes};
    auto pairs = lpm.get_matched_pair_indices();

    ASSERT_EQ(pairs, truth);

}

TEST(CirclePeriodicityMapper, distance_too_big) {

    std::set<std::pair<unsigned int, unsigned int>> truth{
            {1, 5},
            {0, 4},
            {2, 6},
            {3, 7}};

    std::map<unsigned int, std::vector<double>> dof_to_nodes;

    dof_to_nodes[0] = {1.0, 0.0};
    dof_to_nodes[1] = {0.0, 1.0};
    dof_to_nodes[2] = {-1.0, 0.0};
    dof_to_nodes[3] = {0, -1.0};

    double x_offset = 3.5;
    double y_offset = -0.5;
    double eps = 1e-3;
    dof_to_nodes[4] = {1.0+x_offset, 0.0+y_offset};
    dof_to_nodes[5] = {0.0+x_offset, 1.0+y_offset};
    dof_to_nodes[6] = {-1.0+x_offset, 0.0+y_offset+eps};
    dof_to_nodes[7] = {0.0+x_offset, -1.0+y_offset};

    std::vector<unsigned int> a_dofs{1, 0, 2, 3};
    std::vector<unsigned int> b_dofs{4, 5, 6, 7};

    CirclePeriodicityMapper lpm{a_dofs, b_dofs, dof_to_nodes};
    auto pairs = lpm.get_matched_pair_indices();

    ASSERT_EQ(pairs, truth);

}
