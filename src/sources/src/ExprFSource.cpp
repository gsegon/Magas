//
// Created by gordan on 7/13/23.
//

#include "ExprFSource.h"
#include "exprtk.hpp"

ExprFSource::ExprFSource(std::string expression_string) {
    symbol_table.add_variable("x", x);
    symbol_table.add_variable("y", y);
    expression.register_symbol_table(symbol_table);
    parser.compile(expression_string, expression);

}

double ExprFSource::get_value(double x, double y) {
    this->x = x;
    this->y = y;
    return expression.value();
}