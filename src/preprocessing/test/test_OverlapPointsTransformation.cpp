//
// Created by gordan on 3/12/23.
//
#include <gtest/gtest.h>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <fstream>

#include <deal.II/lac/vector.h>
#include <deal.II/base/logstream.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>

#include "OverlapPointsTransformation.h"
#include "vtu11/vtu11.hpp"

using namespace Eigen;
using namespace std;
using namespace dealii;


TEST(OverlapPointsTransformation, basic){


    vector<Vector2f> tocke1{{-1, -1},
                            {1, -1},
                            {1, 1},
                            {-1, 1},
                            {-0.5, 0.5},
                            {0.2, 3}};

    Rotation2D<float> rot2(-M_PI/2/2);
    auto t = Translation<float, 2>(5, 0);

    vector<Vector2f> tocke2(tocke1.size());
    for (int i=0; i < tocke1.size(); i++){
        tocke2[i] = t*rot2*tocke1[i];
    }

    // ----------------------------------
    typedef Point<2> temp_def;
    std::vector<temp_def> a(tocke1.size());
    std::vector<temp_def> b(tocke2.size());

    for (int i=0; i<tocke1.size(); i++){
        a[i] = {tocke1[i][0], tocke1[i][1]};
        b[i] = {tocke2[i][0], tocke2[i][1]};
    }

    std::rotate(b.begin(), b.begin(), b.end());

    std::vector<temp_def> a_trans(tocke1.size());
    OverLapPointsTransformation<temp_def> olpt{a, b};
    olpt.apply_transform(a_trans);

    double sum = 0;
    for (int i=0; i < a_trans.size(); i++){
        std::cout << "b[" << i << "] = " << b[i] << "\t\t\ta_trans[" << i << "] = " << a_trans[i] << std::endl;
        sum += abs(b[i][0]- a_trans[i][0]) + abs(b[i][1]- a_trans[i][1]);
    }
    ASSERT_LT(sum, 1e-6);
}

TEST(OverlapPointsTransformation, motoric_section){

    std::string test_mesh = "/home/gordan/Programs/solver/examples/motoric_section/motoric_section.msh";

    Triangulation<2> triangulation;
    GridIn<2> grid_in;
    grid_in.attach_triangulation(triangulation);
    std::ifstream input_file(test_mesh);
    grid_in.read_msh(input_file);

    DoFHandler<2> dof_handler(triangulation);
    FE_Q<2> fe(1);
    dof_handler.distribute_dofs(fe);

    const IndexSet b_dofs_1 = DoFTools::extract_boundary_dofs(dof_handler, ComponentMask(), {518});
    const IndexSet b_dofs_2 = DoFTools::extract_boundary_dofs(dof_handler, ComponentMask(), {519});

    AssertThrow (b_dofs_1.n_elements() == b_dofs_2.n_elements(), ExcInternalError())

    MappingQ1<2> mapping;
    std::vector<Point<2>> nodes(dof_handler.n_dofs());
    DoFTools::map_dofs_to_support_points(mapping, dof_handler, nodes);

    std::vector<Point<2>> points_a;//(b_dofs_1.n_elements());
    std::vector<Point<2>> b;//(b_dofs_2.n_elements());

    for (auto dof : b_dofs_1){
        points_a.push_back(nodes[dof]);
    }

    for (auto dof : b_dofs_2){
        b.push_back(nodes[dof]);
    }


    ofstream points_a_file;
    points_a_file.open("/home/gordan/Programs/solver/scripts/points_a_file.csv");
    for (auto dof : b_dofs_1){
        points_a_file << nodes[dof][0] << "," << nodes[dof][1] << endl;
    }
    points_a_file.close();

    ofstream points_b_file;
    points_b_file.open("/home/gordan/Programs/solver/scripts/points_b_file.csv");
    for (auto dof : b_dofs_2){
        points_b_file << nodes[dof][0] << "," << nodes[dof][1] << endl;
    }
    points_b_file.close();


    std::vector<Point<2>> a_trans(points_a.size());
    OverLapPointsTransformation<Point<2>> olpt{points_a, b};
    olpt.apply_transform(a_trans);

    double sum = 0;
    for (int i=0; i < a_trans.size(); i++){
        std::cout << "b[" << i << "] = " << b[i] << "\t\t\ta_trans[" << i << "] = " << a_trans[i] << std::endl;
        sum += abs(b[i][0]- a_trans[i][0]) + abs(b[i][1]- a_trans[i][1]);
    }



    ASSERT_LT(sum, 1e-6);

}

TEST(OverlapPointsTransformation, torque_benchmark_kelvin_1){

    std::string test_mesh = "/home/gordan/Programs/solver/examples/torque_benchmark_kelvin_1/torque_benchmark_kelvin_1.msh";

    Triangulation<2> triangulation;
    GridIn<2> grid_in;
    grid_in.attach_triangulation(triangulation);
    std::ifstream input_file(test_mesh);
    grid_in.read_msh(input_file);

    DoFHandler<2> dof_handler(triangulation);
    FE_Q<2> fe(1);
    dof_handler.distribute_dofs(fe);

    const IndexSet b_dofs_1 = DoFTools::extract_boundary_dofs(dof_handler, ComponentMask(), {1});
    const IndexSet b_dofs_2 = DoFTools::extract_boundary_dofs(dof_handler, ComponentMask(), {2});

    AssertThrow (b_dofs_1.n_elements() == b_dofs_2.n_elements(), ExcInternalError())

    MappingQ1<2> mapping;
    std::vector<Point<2>> nodes(dof_handler.n_dofs());
    DoFTools::map_dofs_to_support_points(mapping, dof_handler, nodes);

    std::vector<Point<2>> points_a;//(b_dofs_1.n_elements());
    std::vector<Point<2>> b;//(b_dofs_2.n_elements());

    for (auto dof : b_dofs_1){
        points_a.push_back(nodes[dof]);
    }

    for (auto dof : b_dofs_2){
        b.push_back(nodes[dof]);
    }

    std::vector<Point<2>> a_trans(points_a.size());
    OverLapPointsTransformation<Point<2>> olpt{points_a, b};
    olpt.apply_transform(a_trans);

    double sum = 0;
    for (int i=0; i < a_trans.size(); i++){
        std::cout << "b[" << i << "] = " << b[i] << "\t\t\ta_trans[" << i << "] = " << a_trans[i] << std::endl;
        sum += abs(b[i][0]- a_trans[i][0]) + abs(b[i][1]- a_trans[i][1]);
    }
    ASSERT_LT(sum, 1e-6);


}

TEST(OverlapPointsTransformation, torque_benchmark_kelvin_1_export_points){

    std::string test_mesh = "/home/gordan/Programs/solver/examples/torque_benchmark_kelvin_1/torque_benchmark_kelvin_1.msh";

    Triangulation<2> triangulation;
    GridIn<2> grid_in;
    grid_in.attach_triangulation(triangulation);
    std::ifstream input_file(test_mesh);
    grid_in.read_msh(input_file);

    DoFHandler<2> dof_handler(triangulation);
    FE_Q<2> fe(1);
    dof_handler.distribute_dofs(fe);

    // Points
    std::vector<double> points_a;
    std::vector<double> points_b;

    const IndexSet b_dofs_1 = DoFTools::extract_boundary_dofs(dof_handler, ComponentMask(), {1});
    const IndexSet b_dofs_2 = DoFTools::extract_boundary_dofs(dof_handler, ComponentMask(), {2});

    AssertThrow (b_dofs_1.n_elements() == b_dofs_2.n_elements(), ExcInternalError())

    MappingQ1<2> mapping;
    std::vector<Point<2>> nodes(dof_handler.n_dofs());
    DoFTools::map_dofs_to_support_points(mapping, dof_handler, nodes);

    std::vector<Point<2>> a, b;

    for (auto dof : b_dofs_1){
        points_a.push_back(nodes[dof][0]);
        points_a.push_back(nodes[dof][1]);
        points_a.push_back(0.0);
    }

    for (auto dof : b_dofs_2){
        points_b.push_back(nodes[dof][0]);
        points_b.push_back(nodes[dof][1]);
        points_b.push_back(0.0);
    }

    ofstream points_a_file;
    points_a_file.open("/home/gordan/Programs/solver/scripts/points_a_file.csv");
    for (auto dof : b_dofs_1){
        points_a_file << nodes[dof][0] << "," << nodes[dof][1] << endl;
    }
    points_a_file.close();

    ofstream points_b_file;
    points_b_file.open("/home/gordan/Programs/solver/scripts/points_b_file.csv");
    for (auto dof : b_dofs_2){
        points_b_file << nodes[dof][0] << "," << nodes[dof][1] << endl;
    }
    points_b_file.close();

}

