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

#include <gtest/gtest.h>
#include <cstdio>
#include <fstream>

#include "strtk.hpp"

#include "InterpolatedNuCurve.h"

TEST(InterpolatedBHCurve, basic) {

    std::vector<double> b{0, 0.2003, 0.3204, 0.40045, 0.50055, 0.5606,
                          0.7908, 0.931, 1.1014,
                          1.2016, 1.302, 1.4028,
                          1.524, 1.626, 1.698,
                          1.73, 1.87, 1.99,
                          2.04, 2.07, 2.095,
                          2.2, 2.4};

    std::vector<double> h{10, 238.7, 318.3,
                          358.1, 437.7, 477.5,
                          636.6, 795.8, 1114.1,
                          1273.2, 1591.5, 2228.2,
                          3183.1, 4774.6, 6366.2,
                          7957.7, 15915.5, 31831,
                          47746.5, 63663, 79577.5,
                          159155, 318310};

    InterpolatedNuCurve ibh(b, h);

//    std::cout << "nu(0.15): " << ibh.get_nu(0.15) << std::endl;
//    std::cout << "nu_prime(0.15): " << ibh.get_nu_prime(0.15) << std::endl;

    auto fbh = fopen("bh.dat", "w");
    auto fnu = fopen("nu.dat", "w");
    for (double bi=b[0]; bi < b[b.size()-1]; bi += (b[b.size()-1]-b[0])/100){
        double nui = ibh.get_nu(bi);
        double hi = nui*bi;
        fprintf(fbh, "%g %g\n", hi, bi);
        fprintf(fnu, "%g %g\n", bi, nui);
    }
    fclose(fbh);
    fclose(fnu);


}

TEST(InterpolatedBHCurve, from_csv) {

    std::vector<double> b;
    std::vector<double> h;

    std::string filename{"../../../src/material/test/M-45.csv"};
    strtk::token_grid::options options;
    options.column_delimiters =",";

    strtk::token_grid bh_grid(filename, options);
    for (std::size_t r = 0; r < bh_grid.row_count(); ++r){
        strtk::token_grid::row_type row = bh_grid.row(r);

        if (r==0){
        }
        else{
            b.push_back(std::stod(row.get<std::string>(0)));
            h.push_back(std::stod(row.get<std::string>(1)));
        }
    }

    InterpolatedNuCurve ibh(b, h);
    std::cout << "nu(0): " << ibh.get_nu(0) << std::endl;
    std::cout << "nu_prime(0): " << ibh.get_nu_prime(0.0) << std::endl;

    std::cout << "nu(b[b.size()]+1): " << ibh.get_nu(b[b.size()-1]+1) << std::endl;
    std::cout << "nu_prime(b[b.size()-1]+1): " << ibh.get_nu_prime(b[b.size()-1]+1) << std::endl;

    auto fbh = fopen("bh.dat", "w");
    auto fnu = fopen("nu.dat", "w");
    auto fbce = fopen("bce.dat", "w");
    auto fbe = fopen("be.dat", "w");
    for (double bi=b[0]; bi < 1.5*b[b.size()-1]; bi += (1.5*b[b.size()-1]-b[0])/500.0){
        double nui = ibh.get_nu(bi);
        double hi = nui*bi;
        double ce = ibh.get_coenergy(bi);
        double e = ibh.get_energy(bi);
        fprintf(fbh, "%g %g\n", hi, bi);
        fprintf(fnu, "%g %g\n", bi, nui);
        fprintf(fbce, "%g %g\n", bi, ce);
        fprintf(fbe, "%g %g\n", bi, e);
    }
    fclose(fbh);
    fclose(fnu);
    fclose(fbce);
    fclose(fbe);

}
