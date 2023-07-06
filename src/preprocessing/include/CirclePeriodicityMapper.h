//
// Created by gordan on 6/15/23.
//

#ifndef SOLVER_CIRCLEPERIODICITYMAPPER_H
#define SOLVER_CIRCLEPERIODICITYMAPPER_H

#include <vector>
#include <limits>
#include <cmath>
#include <map>
#include <set>
#include <cassert>
#include <algorithm>
#include <cmath>
#include "IPeriodicityMapper.h"


using namespace std;

class CirclePeriodicityMapper : public IPeriodicityMapper{

public:
    CirclePeriodicityMapper(vector<unsigned int>,
                            vector<unsigned int>,
                            map<unsigned int, vector<double>>);

    set<pair<unsigned int, unsigned int>> get_matched_pair_indices() override;
    map<unsigned int, vector<double>> reference_to_center(map<unsigned int, vector<double>>);
    double get_weigth();
    void set_weigth(double);
    double error_overlap();

};


#endif //SOLVER_CIRCLEPERIODICITYMAPPER_H
