//
// Created by gordan on 7/3/23.
//

#include "processors/ScalarPostprocessorFactory.h"
#include "processors/ExpressionScalarPostprocessor.h"

template class ScalarPostprocessorFactory<2>;

template<int dim>
ScalarPostprocessorFactory<dim>::ScalarPostprocessorFactory(const std::unordered_map<int, double>& nu_map, std::unordered_map<int, std::variant<double, std::pair<double, double>>>& f_map){
    this->nu_map_ptr = &nu_map;
    this->f_map_ptr = &f_map;
}

template<int dim>
ScalarPostprocessor<dim> *ScalarPostprocessorFactory<dim>::create(std::string input_string) {
    return new ExpressionScalarPostprocessor<dim>(input_string, *nu_map_ptr, *f_map_ptr);
}