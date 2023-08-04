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
// Created by gordan on 7/3/23.
//

#ifndef MAGAS_SCALARPOSTPROCESSORFACTORY_H
#define MAGAS_SCALARPOSTPROCESSORFACTORY_H

#include <unordered_map>
#include <variant>

#include "ScalarPostprocessor.h"
#include "NuCurve.h"
#include "FSource.h"

template<int dim>
class ScalarPostprocessorFactory {

public:
    ScalarPostprocessorFactory(const std::unordered_map<int, NuCurve*>&, std::unordered_map<int, std::variant<FSource*, std::pair<double, double>>>&);
    ScalarPostprocessor<dim>* create(std::string input_string);

protected:
    const std::unordered_map<int, NuCurve*>* nu_map_ptr = nullptr;
    std::unordered_map<int, std::variant<FSource*, std::pair<double, double>>>* f_map_ptr = nullptr;

};


#endif //MAGAS_SCALARPOSTPROCESSORFACTORY_H
