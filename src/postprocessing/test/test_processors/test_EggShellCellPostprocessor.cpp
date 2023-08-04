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

//
// Created by gordan on 3/12/23.
//
#include <gtest/gtest.h>
#include <tuple>
#include <unordered_map>
#include <string>
#include <any>
#include <fstream>
#include <filesystem>

#include "LinearSolver.h"
#include "NewtonSolver.h"
#include "processors/EggShellCellPostprocessor.h"
#include "processors/MatIDPostprocessor.h"
#include "NuCurve.h"
#include "LinearNuCurve.h"
#include "export/ExportVtu.h"
#include "ConstFSource.h"

TEST(TestEggShellCellPostprocessor, torque_benchmark_kelvin_1){

    std::filesystem::path home = std::getenv("HOME");
    std::filesystem::path test_mesh = "../../../examples/torque_benchmark_kelvin_1/torque_benchmark_kelvin_1.msh";


    std::unordered_map<int, NuCurve*> nu_map{{3, new LinearNuCurve{795774.715025}},
                                             {4, new LinearNuCurve{795774.715025}},
                                             {5, new LinearNuCurve{795774.715025}},
                                             {6, new LinearNuCurve{795774.715025}},
                                             {8, new LinearNuCurve{795774.715025}}};

    std::unordered_map<int, std::variant<FSource*, std::pair<double, double>>> f_map{{3, std::pair<double, double>{0.0, 1000000.0}},
                                                                                   {4, new ConstFSource{0}},
                                                                                   {5, new ConstFSource{0}},
                                                                                   {6, new ConstFSource{0}},
                                                                                   {8, std::pair<double, double>{-1591549.43091895, 0.0}}
    };
    std::unordered_map<int, double> dc_map{{7, 0}};

    std::unordered_map<std::string, std::vector<unsigned int>> per_map{{"periodic-circle", {1, 2}}};

    LinearSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.set_nu_map(nu_map);
    solver.set_f_map(f_map);
    solver.set_dc_map(dc_map);
    solver.set_per_map(per_map);
    solver.setup_system();
    solver.assemble_system();
    solver.solve();

    MatIDPostprocessor<2> mat_id_postp{};
    EggShellCellPostprocessor<2> neighbor_postp{3, 4};

    ExportVtu<2> export_vtu(solver.get_triangulation(), solver.get_rhs(), solver.get_solution(), solver.get_fe());
    export_vtu.attach_postprocessor(&neighbor_postp, "Eggshell");
    export_vtu.attach_postprocessor(&mat_id_postp, "mat_id");

    export_vtu.write("neighbors_torque_benchmark_kelvin_1");

}


