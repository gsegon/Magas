//
// Created by gordan on 7/5/23.
//

#include "../include/InterpolatedNuCurve.h"
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <iostream>
#include <cassert>

InterpolatedNuCurve::InterpolatedNuCurve(std::vector<double> b, std::vector<double> h) {
   this->b = b;
   this->h = h;

   if (this->b.size() != this->h.size()){
       throw std::runtime_error("B entry not the same size as H entry.");
   }

   this->acc = gsl_interp_accel_alloc();
   this->spline = gsl_spline_alloc(gsl_interp_cspline, b.size());

   gsl_spline_init(spline, &b[0], &h[0], b.size());
}

InterpolatedNuCurve::~InterpolatedNuCurve() {
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
}

double InterpolatedNuCurve::get_nu(double b) {
    if (b!=0){
        double h = gsl_spline_eval(spline, b, acc);
        return h/b;
    }
    else{
        return InterpolatedNuCurve::get_nu(1e-3);
    }

}

double InterpolatedNuCurve::get_nu_prime(double b) {
    if (b!=0){
        double h = gsl_spline_eval(spline, b, acc);
        double h_prime = gsl_spline_eval_deriv(spline, b, acc);
        return h_prime/b - h/(b*b);
    }
    else{
        return InterpolatedNuCurve::get_nu_prime(1e-3);
    }
}
