// Magas - Magnetostatic Analysis Suite
// Copyright (C) 2023  Gordan Segon
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
// USA

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
    bool is_full_circle();
    void print_map();
    dof get_other(dof);
    bool is_in(dof);

private:
    void sort();

    map<dof, vector<double>> dof_to_node;
    vector<dof> dofs;
    int offset;
    bool full_circle;

};


#endif //MAGAS_SLIDINGROTATION_H
