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

#ifndef MAGAS_INTERPOLATEDNUCURVE_H
#define MAGAS_INTERPOLATEDNUCURVE_H

#include <vector>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "NuCurve.h"

class InterpolatedNuCurve : public NuCurve {

public:
    InterpolatedNuCurve(std::vector<double>, std::vector<double>);
    ~InterpolatedNuCurve();
    double get_nu(double) override;
    double get_nu_prime(double) override;
    double get_coenergy(double);
    double get_energy(double);

private:
    std::vector<double> b;
    std::vector<double> h;

    gsl_interp_accel *acc = nullptr;
    gsl_spline *spline = nullptr;

};


#endif //MAGAS_INTERPOLATEDNUCURVE_H
