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
