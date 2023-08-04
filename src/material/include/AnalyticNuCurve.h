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

#ifndef MAGAS_ANALYTICNUCURVE_H
#define MAGAS_ANALYTICNUCURVE_H

#include "NuCurve.h"

class AnalyticNuCurve : public NuCurve {

public:
    AnalyticNuCurve();
    double get_nu(double) override;
    double get_nu_prime(double) override;
    double get_energy(double) override;
    double get_coenergy(double) override;

private:
    double nu_0;
    double alpha;
    double beta;
    double theta;
    double b_sat;
};


#endif //MAGAS_ANALYTICNUCURVE_H
