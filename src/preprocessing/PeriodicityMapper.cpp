//
// Created by gordan on 6/15/23.
//
#include <vector>
#include <limits>
#include <cmath>
#include <map>
#include <cassert>

#include "PeriodicityMapper.h"

template class PeriodicityMapper<std::tuple<double, double>>;

template<typename T>
PeriodicityMapper<T>::PeriodicityMapper(const std::vector<T>& firsts, const std::vector<T>& seconds){

        assert(firsts.size() == seconds.size() && "firsts size not equal to seconds.");

        this->firsts = firsts;
        this->seconds = seconds;

        for (auto val : firsts){
            x_max = (std::get<0>(val) > x_max) ? std::get<0>(val) : x_max;
            x_min = (std::get<0>(val) < x_min) ? std::get<0>(val) : x_min;
            y_max = (std::get<1>(val) > y_max) ? std::get<1>(val) : y_max;
            y_min = (std::get<1>(val) < y_min) ? std::get<1>(val) : y_min;
            s_max = (x_max > y_max) ? x_max : y_max;
            s_min = (x_min > y_min) ? x_min : y_min;
        }

        s_span = s_max - s_min;
        n_x = double(sqrt(std::numeric_limits<unsigned int>::max()));
        n_y = double(sqrt(std::numeric_limits<unsigned int>::max()));

}

template<typename T>
unsigned int PeriodicityMapper<T>::hash(T point){
    auto x = std::get<0>(point);
    auto y = std::get<1>(point);
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
