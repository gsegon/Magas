//
// Created by gordan on 3/12/23.
//
#include <gtest/gtest.h>

#include "PeriodicityMapper.h"

TEST(PeriodicityMapper, basic){

    std::vector<std::vector<double>> firsts{{0, 0},
                                                   {0, 1},
                                                   {0.3, 0.2},
                                                   {0.5, 0.5},
                                                   {1, 1}};

    std::vector<std::vector<double>> seconds{{0, 1},
                                                    {0.3, 0.2},
                                                    {0, 0},
                                                    {0.5, 0.5},
                                                    {1, 1}};
    PeriodicityMapper<std::vector<double>> pm{firsts, seconds};

    pm.map_points();
    auto matched_pairs = pm.get_matched_pair_indices();

    std::cout << "Matched pairs: " << std::endl;
    for (auto pair : matched_pairs)
        std::cout << "(" << pair.first << ", " << pair.second << ")" << std::endl;

}

