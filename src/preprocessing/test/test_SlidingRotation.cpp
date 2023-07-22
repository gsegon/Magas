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