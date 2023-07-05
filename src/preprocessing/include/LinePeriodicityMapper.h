//
// Created by gordan on 6/15/23.
//

#ifndef SOLVER_LINEPERIODICITYMAPPER_H
#define SOLVER_LINEPERIODICITYMAPPER_H

#include <vector>
#include <limits>
#include <cmath>
#include <map>
#include <cassert>

#include "IPeriodicityMapper.h"

using namespace std;

class LinePeriodicityMapper : public IPeriodicityMapper{

public:
    LinePeriodicityMapper(vector<unsigned int>,
                            vector<unsigned int>,
                            map<unsigned int, vector<double>>);

    set<pair<unsigned int, unsigned int>> get_matched_pair_indices() override;
    double get_weigth() override;
    void set_weigth(double);

};


#endif //SOLVER_LINEPERIODICITYMAPPER_H
