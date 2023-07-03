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

    std::vector<std::string> patterns{arkkio, b_abs, lorentz, flux_linkage};
}

#endif //MAGAS_REGEXPATTERNS_H
