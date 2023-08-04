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

#ifndef MAGAS_REGEXPATTERNS_H
#define MAGAS_REGEXPATTERNS_H

#include <string>
#include <vector>

namespace RegexPatterns{

    std::string arkkio{R"((Arkkio)\((.*),\s*\((.*),(.*)\)\))"};
    std::string b_abs{R"((Babs)\((.*),\s*(.*)\))"};
    std::string lorentz{R"((LorentzForce)\((.*),\s*(.*)\))"};
    std::string flux_linkage{R"((Psi)\((.*)\))"};
    std::string torque_eggshell{R"((Torque)\((.*),\s*(.*)\))"};
    std::string force_eggshell{R"((Force)\((.*),\s*(.*),\s*(.*)\))"};

    std::vector<std::string> patterns{arkkio, b_abs, lorentz, flux_linkage, torque_eggshell, force_eggshell};
}

#endif //MAGAS_REGEXPATTERNS_H
