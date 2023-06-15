//
// Created by gordan on 6/15/23.
//

#ifndef SOLVER_PERIODICITYMAPPER_H
#define SOLVER_PERIODICITYMAPPER_H

#include <vector>
#include <limits>
#include <cmath>
#include <map>

template<typename T>
class PeriodicityMapper{

public:
    PeriodicityMapper(const std::vector<T>& firsts, const std::vector<T>& seconds);
    unsigned int hash(T point);
    void map_points();
    std::vector<std::pair<unsigned int, unsigned int>> get_matched_pair_indices();


private:
    std::vector<T> firsts;
    std::vector<T> seconds;
    std::map<unsigned int, std::pair<int, int>> periodicity_map;
    double x_max = std::numeric_limits<double>::min();
    double x_min = std::numeric_limits<double>::max();
    double y_max = std::numeric_limits<double>::min();
    double y_min = std::numeric_limits<double>::max();
    double s_max = std::numeric_limits<double>::min();
    double s_min = std::numeric_limits<double>::max();

    double s_span;
    double n_x;
    double n_y;
};


#endif //SOLVER_PERIODICITYMAPPER_H
