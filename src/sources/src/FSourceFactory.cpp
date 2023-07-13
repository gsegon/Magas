//
// Created by gordan on 7/13/23.
//

#include "FSourceFactory.h"
#include "ConstFSource.h"
#include "ExprFSource.h"

FSource* FSourceFactory::create(double val) {
    return new ConstFSource{val};
}

FSource* FSourceFactory::create(std::string expr_string) {
    return new ExprFSource{expr_string};
}