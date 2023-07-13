//
// Created by gordan on 7/13/23.
//

#include "../include/ConstFSource.h"

ConstFSource::ConstFSource(double value) {
    this->value = value;
}

double ConstFSource::get_value(double x, double y) {
    return value;
}