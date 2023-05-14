//
// Created by gordan on 3/12/23.
//
#include <gtest/gtest.h>
#include <tuple>
#include <unordered_map>
#include <string>
#include <any>
#include <fstream>

#include "LinearSolver.h"
#include "PicardSolver.h"
#include "PostMagneticFluxDensity.h"
#include "vtu11/vtu11.hpp"

TEST(TestVtu11, triangulation_to_vtu){

    std::string test_mesh = "/home/gordan/Programs/solver/test/test_data/test_EI_core/EI_core.msh";

    Triangulation<2> triangulation;
    GridIn<2> grid_in;
    grid_in.attach_triangulation(triangulation);
    std::ifstream input_file(test_mesh);
    grid_in.read_msh(input_file);

    // Some cell data
    std::vector<double> mat_ids;

    // Points
    std::vector<double> points;
    for(const auto& vertex : triangulation.get_vertices()){
        points.push_back(vertex[0]); points.push_back(vertex[1]); points.push_back(0.0);
    }

    std::vector<vtu11::VtkIndexType> connectivity;
    for(const auto& cell: triangulation.active_cell_iterators()){
        for (auto local_vertex_index: cell->vertex_indices()){
            connectivity.push_back(cell->vertex_index(local_vertex_index));
        }
        std::cout << cell->material_id() << std::endl;
        mat_ids.push_back(cell->material_id());
    }

    std::vector<vtu11::VtkIndexType> offsets(triangulation.n_active_cells());
    long n = 0;
    std::generate(offsets.begin(), offsets.end(), [&n]{ return n+=4; });

    std::vector<vtu11::VtkCellType> types(triangulation.n_active_cells(), 8);

    vtu11::Vtu11UnstructuredMesh mesh { points, connectivity, offsets, types };

    // Create tuples with (name, association, number of components) for each data set
    std::vector<vtu11::DataSetInfo> dataSetInfo{
            { "Material ID", vtu11::DataSetType::CellData, 1},
    };



    vtu11::writeVtu("test_2.vtu", mesh, dataSetInfo, {mat_ids}, "Ascii");


}


