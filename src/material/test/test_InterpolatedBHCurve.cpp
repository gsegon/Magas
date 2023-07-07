//
// Created by gordan on 3/12/23.
//
#include <gtest/gtest.h>
#include <cstdio>

#include "../include/InterpolatedBHCurve.h"

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

    InterpolatedBHCurve ibh(b, h);

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
