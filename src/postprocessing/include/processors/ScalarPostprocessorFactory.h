//
// Created by gordan on 7/3/23.
//

#ifndef MAGAS_SCALARPOSTPROCESSORFACTORY_H
#define MAGAS_SCALARPOSTPROCESSORFACTORY_H

#include <unordered_map>
#include <variant>

#include "ScalarPostprocessor.h"

template<int dim>
class ScalarPostprocessorFactory {

public:
    ScalarPostprocessorFactory(const std::unordered_map<int, double>&, std::unordered_map<int, std::variant<double, std::pair<double, double>>>&);
    ScalarPostprocessor<dim>* create(std::string input_string);

protected:
    const std::unordered_map<int, double>* nu_map_ptr = nullptr;
    std::unordered_map<int, std::variant<double, std::pair<double, double>>>* f_map_ptr = nullptr;

};


#endif //MAGAS_SCALARPOSTPROCESSORFACTORY_H
