//
// Created by gordan on 7/11/23.
//

#ifndef MAGAS_SLIDINGROTATION_H
#define MAGAS_SLIDINGROTATION_H

#include <vector>
#include <map>

using namespace std;

typedef unsigned int dof;

class SlidingRotation {

public:
    SlidingRotation(vector<dof> dofs, map<dof, vector<double>> dof_to_node, int offset);
    dof get_mapped(dof);
    vector<dof> get_dofs();

private:
    void sort();

    map<dof, vector<double>> dof_to_node;
    vector<dof> dofs;
    int offset;

};


#endif //MAGAS_SLIDINGROTATION_H
