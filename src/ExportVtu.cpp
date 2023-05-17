//
// Created by gordan on 5/15/23.
//

#include <fstream>
#include "ExportVtu.h"
#include "vtu11/vtu11.hpp"


#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>



template class ExportVtu<2>;

template <int dim>
ExportVtu<dim>::ExportVtu(const Triangulation<dim>& triangulation) {
    triangulation_ptr = &triangulation;
}


template <int dim>
ExportVtu<dim>::ExportVtu(const Triangulation<dim> &triangulation, const Vector<double>& rhs,
                          const Vector<double>& solution, const FE_Q<dim>& fe) {
    triangulation_ptr = &triangulation;
    rhs_ptr = &rhs;
    solution_ptr = &solution;
    fe_ptr = &fe;
}


template <int dim>
void ExportVtu<dim>::write(const std::string& filename){

    DoFHandler<dim> dof_handler(*triangulation_ptr);
    dof_handler.distribute_dofs(*fe_ptr);

    // Cell data
    std::vector<double> mat_ids;

    // Point data
    std::vector<double> rhs_export;
    std::vector<double> solution_export;

    // Points
    std::vector<double> points;

    int point_index = 0;
    std::vector<vtu11::VtkIndexType> connectivity;
    for (auto cell : dof_handler.active_cell_iterators()){
        for (unsigned int vertex_index = 0; vertex_index < GeometryInfo<2>::vertices_per_cell; ++vertex_index){
            for (unsigned int component=0; component<fe_ptr->n_components(); ++component){

                const unsigned int dof_index = cell->vertex_dof_index(vertex_index,component);

                // Add vertex to points
                auto vertex = cell->vertex(vertex_index);
                points.push_back(vertex[0]); points.push_back(vertex[1]); points.push_back(0.0);

                // Solution and rhs
                solution_export.push_back((*solution_ptr)[dof_index]);
                rhs_export.push_back((*rhs_ptr)[dof_index]);

            }
            connectivity.push_back(point_index);
            connectivity.push_back(point_index+1);
            connectivity.push_back(point_index+3);
            connectivity.push_back(point_index+2);
            point_index+=4;

        }
        // Mat id data
        mat_ids.push_back(cell->material_id());
    }

    std::vector<vtu11::VtkIndexType> offsets(triangulation_ptr->n_active_cells());
    long n = 0;
    std::generate(offsets.begin(), offsets.end(), [&n]{ return n+=4; });

    std::vector<vtu11::VtkCellType> types(triangulation_ptr->n_active_cells(), 9);
    vtu11::Vtu11UnstructuredMesh mesh { points, connectivity, offsets, types };


    std::vector<std::tuple<double, double>> B;
    std::vector<double> Bx;
    std::vector<double> By;

    QGauss<dim> quadrature_formula(fe_ptr->degree +1);

    FEValues<dim> fe_values(*fe_ptr, quadrature_formula, update_values | update_gradients | update_quadrature_points |
                                                         update_JxW_values);
    std::vector<Tensor<1, dim>> solution_gradients(quadrature_formula.size());
    for (auto cell : dof_handler.active_cell_iterators()){
        fe_values.reinit(cell);
        fe_values.get_function_gradients(*solution_ptr, solution_gradients);

        Bx.push_back(solution_gradients[0][1]);
        By.push_back(-solution_gradients[0][0]);

        B.push_back(std::tuple<double, double>{solution_gradients[0][1], -solution_gradients[0][0]});
    }



    // Create tuples with (name, association, number of components) for each data set
    std::vector<vtu11::DataSetInfo> dataSetInfo{
            { "Material ID", vtu11::DataSetType::CellData, 1},
            { "rhs", vtu11::DataSetType::PointData, 1},
            { "solution", vtu11::DataSetType::PointData, 1},
            { "Bx", vtu11::DataSetType::CellData, 1},
            { "By", vtu11::DataSetType::CellData, 1},
            { "B", vtu11::DataSetType::CellData, 2},
    };

    vtu11::writeVtu(filename+".vtu", mesh, dataSetInfo, {mat_ids, rhs_export, solution_export, Bx, By, Bx},"Ascii");

}

