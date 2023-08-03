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
// Created by gordan on 7/13/23.
//

#ifndef MAGAS_EXPRFSOURCE_H
#define MAGAS_EXPRFSOURCE_H

#include <string>
#include "exprtk.hpp"

#include "FSource.h"


typedef exprtk::symbol_table<double> symbol_table_t;
typedef exprtk::expression<double>   expression_t;
typedef exprtk::parser<double>       parser_t;

class ExprFSource : public FSource {

public:
    ExprFSource(std::string expression_string);
    double get_value(double x, double y);

private:
    std::string expression_string;
    symbol_table_t symbol_table;
    expression_t expression;
    parser_t parser;
    double x;
    double y;

};


#endif //MAGAS_EXPRFSOURCE_H
