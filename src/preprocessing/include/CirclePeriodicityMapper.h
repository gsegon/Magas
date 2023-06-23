//
// Created by gordan on 6/15/23.
//

#ifndef SOLVER_CIRCLEPERIODICITYMAPPER_H
#define SOLVER_CIRCLEPERIODICITYMAPPER_H

//#include <vector>
//#include <limits>
//#include <cmath>
//#include <map>
//#include "IPeriodicityMapper.h"
//
//template<typename T>
//class CirclePeriodicityMapper : public IPeriodicityMapper{
//
//public:
//    CirclePeriodicityMapper(std::vector<unsigned int>,
//                            std::vector<unsigned int>,
//                            std::vector<std::vector<double>>);
//    unsigned int hash(T point);
//    void map_points();
//    std::set<std::pair<unsigned int, unsigned int>> get_matched_pair_indices();
//
//
//private:
//    std::vector<T> firsts;
//    std::vector<T> seconds;
//    std::map<unsigned int, std::pair<int, int>> periodicity_map;
//    double x_max = std::numeric_limits<double>::min();
//    double x_min = std::numeric_limits<double>::max();
//    double y_max = std::numeric_limits<double>::min();
//    double y_min = std::numeric_limits<double>::max();
//    double s_max = std::numeric_limits<double>::min();
//    double s_min = std::numeric_limits<double>::max();
//
//    double s_span;
//    double n_x;
//    double n_y;
//};
//
//template<typename T>
//CirclePeriodicityMapper<T>::CirclePeriodicityMapper(std::vector<unsigned int> a_dofs,
//                                                    std::vector<unsigned int> b_dofs,
//                                                    std::vector<std::vector<double>> nodes){
//
//}
//
//template<typename T>
//unsigned int CirclePeriodicityMapper<T>::hash(T point){
//    auto x = point[0];
//    auto y = point[1];
//    // TODO: Reimplement with ints only.
//    return (unsigned int)(  ((x-s_min)/s_span)*n_x + (y-s_min)/(s_span/n_y) *(n_x-1));
//}
//
//template<typename T>
//void CirclePeriodicityMapper<T>::map_points(){
//    for (int index=0; index < firsts.size(); index++){
//        if (!periodicity_map.count(hash(firsts[index])))
//            periodicity_map[hash(firsts[index])] = {};
//
//        if (!periodicity_map.count(hash(seconds[index])))
//            periodicity_map[hash(seconds[index])] = {};
//
//        periodicity_map[hash(firsts[index])].first = index;
//        periodicity_map[hash(seconds[index])].second = index;
//    }
//    assert(periodicity_map.size() == firsts.size() && "periodicity map size not equal to firsts size.");
//}
//
//template<typename T>
//std::set<std::pair<unsigned int, unsigned int>> CirclePeriodicityMapper<T>::get_matched_pair_indices(){
//    std::vector<std::pair<unsigned int, unsigned int>> matched_pairs;
//    matched_pairs.reserve(periodicity_map.size());
//    for (auto& [key, value] : periodicity_map){
//        matched_pairs.emplace_back(value);
//    }
//    return matched_pairs;
//}


#endif //SOLVER_CIRCLEPERIODICITYMAPPER_H
