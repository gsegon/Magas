//
// Created by gordan on 7/3/23.
//
#include <regex>

#include "processors/ScalarPostprocessorFactory.h"
#include "processors/ExpressionScalarPostprocessor.h"
#include "processors/RegexPatterns.h"

#include "processors/ArkkioScalarPostprocessor.h"
#include "processors/PointBabsScalarPostprocessor.h"
#include "processors/FluxLinkageScalarPostprocessor.h"
#include "processors/TorqueEggShellScalarPostprocessor.h"
#include "processors/ForceEggShellScalarPostprocessor.h"
#include "NuCurve.h"

template class ScalarPostprocessorFactory<2>;

template<int dim>
ScalarPostprocessorFactory<dim>::ScalarPostprocessorFactory(const std::unordered_map<int, NuCurve*>& nu_map, std::unordered_map<int, std::variant<FSource*, std::pair<double, double>>>& f_map){
    this->nu_map_ptr = &nu_map;
    this->f_map_ptr = &f_map;
}

template<int dim>
ScalarPostprocessor<dim> *ScalarPostprocessorFactory<dim>::create(std::string input_string) {

    for (auto pattern : RegexPatterns::patterns) {
        std::smatch sm;
        if (std::regex_match(input_string, sm, std::regex(pattern))) {
            if (sm[1] == "Arkkio") return new ArkkioScalarPostprocessor<dim>(std::stoi(sm[2]), *nu_map_ptr);
            if (sm[1] == "Babs") return new PointBabsScalarPostprocessor<dim>({std::stod(sm[2]), std::stod(sm[3])});
            if (sm[1] == "Psi") return new FluxLinkageScalarPostprocessor<dim>(std::stoi(sm[2]), *f_map_ptr);
            if (sm[1] == "Torque") return new TorqueEggShellScalarPostprocessor<dim>({std::stoi(sm[2]), std::stoi(sm[3])});
            if (sm[1] == "Force") return new ForceEggShellScalarPostprocessor<dim>({std::stoi(sm[2]), std::stoi(sm[3]), sm[4]});
        }
    }
    return new ExpressionScalarPostprocessor<dim>(input_string, *nu_map_ptr, *f_map_ptr);
}