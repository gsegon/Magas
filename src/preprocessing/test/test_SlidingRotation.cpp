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