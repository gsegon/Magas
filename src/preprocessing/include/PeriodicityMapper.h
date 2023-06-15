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

template<typename T>
PeriodicityMapper<T>::PeriodicityMapper(const std::vector<T>& firsts, const std::vector<T>& seconds){

    assert(firsts.size() == seconds.size() && "firsts size not equal to seconds.");

    this->firsts = firsts;
    this->seconds = seconds;

    for (auto val : firsts){
        x_max = (val[0] > x_max) ? val[0] : x_max;
        x_min = (val[0] < x_min) ? val[0] : x_min;
        y_max = (val[1] > y_max) ? val[1] : y_max;
        y_min = (val[1] < y_min) ? val[1] : y_min;
        s_max = (x_max > y_max) ? x_max : y_max;
        s_min = (x_min > y_min) ? x_min : y_min;
    }

    s_span = s_max - s_min;
    n_x = double(sqrt(std::numeric_limits<unsigned int>::max()));
    n_y = double(sqrt(std::numeric_limits<unsigned int>::max()));

}

template<typename T>
unsigned int PeriodicityMapper<T>::hash(T point){
    auto x = point[0];
    auto y = point[1];
    // TODO: Reimplement with ints only.
    return (unsigned int)(  ((x-s_min)/s_span)*n_x + (y-s_min)/(s_span/n_y) *(n_x-1));
}

template<typename T>
void PeriodicityMapper<T>::map_points(){
    for (int index=0; index < firsts.size(); index++){
        if (!periodicity_map.count(hash(firsts[index])))
            periodicity_map[hash(firsts[index])] = {};

        if (!periodicity_map.count(hash(seconds[index])))
            periodicity_map[hash(seconds[index])] = {};

        periodicity_map[hash(firsts[index])].first = index;
        periodicity_map[hash(seconds[index])].second = index;
    }
    assert(periodicity_map.size() == firsts.size() && "periodicity map size not equal to firsts size.");
}

template<typename T>
std::vector<std::pair<unsigned int, unsigned int>> PeriodicityMapper<T>::get_matched_pair_indices(){
    std::vector<std::pair<unsigned int, unsigned int>> matched_pairs;
    matched_pairs.reserve(periodicity_map.size());
    for (auto& [key, value] : periodicity_map){
        matched_pairs.emplace_back(value);
    }
    return matched_pairs;
}


#endif //SOLVER_PERIODICITYMAPPER_H
