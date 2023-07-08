//
// Created by gordan on 7/5/23.
//

#include "InterpolatedNuCurve.h"
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

    if (b > this->b[this->b.size()-1]){
        double diff_b = b-this->b[this->b.size()-1];
        double diff_h = 795774.715025*diff_b;
        double h_new = this->h[this->h.size()-1] + diff_h;
        return h_new/b;
    }
    else if (b!=0){
        double h = gsl_spline_eval(spline, b, acc);
        return h/b;
    }
    else{
        return InterpolatedNuCurve::get_nu(1e-3);
    }
}

double InterpolatedNuCurve::get_nu_prime(double b) {

    if (b > this->b[this->b.size()-1]){
        return 0;
    }
    else if (b!=0){
        double h = gsl_spline_eval(spline, b, acc);
        double h_prime = gsl_spline_eval_deriv(spline, b, acc);
        return h_prime/b - h/(b*b);
    }
    else if(b==0){
        return 0;
    }
}

double InterpolatedNuCurve::get_coenergy(double b) {
    double coenergy;
    if (b==0){
        coenergy = 0;
    }
    else if (b > this->b[this->b.size()-1]){
        double e_nonlinear_part = gsl_spline_eval_integ(spline,0, this->b[this->b.size()-1], acc);
        double b_interp_max = this->b[this->b.size()-1];
        double h_interp_max = this->h[this->h.size()-1];
        double diff_b = b-b_interp_max;
        double diff_h = 795774.715025*diff_b;
        double e_linear_part = diff_b*(h_interp_max+diff_h/2.0);
        coenergy = e_nonlinear_part + e_linear_part;
    }
    else{
        coenergy = gsl_spline_eval_integ(spline, 0, b, acc);
    }

    return coenergy;
}

double InterpolatedNuCurve::get_energy(double b) {
    double energy;

    if (b==0){
        energy = 0;
    }
    else if (b > this->b[this->b.size()-1]){
        double diff_b = b-this->b[this->b.size()-1];
        double diff_h = 795774.715025*diff_b;
        double h_new = this->h[this->h.size()-1] + diff_h;
        energy = h_new*b - get_coenergy(b);
    }
    else{
        double h_new = get_nu(b)*b;
        energy = h_new*b - get_coenergy(b);
    }

    return energy;
}