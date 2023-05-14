//
// Created by gordan on 5/14/23.
//

#ifndef SOLVER_POSTMAGNETICFLUXDENSITY_H
#define SOLVER_POSTMAGNETICFLUXDENSITY_H

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/patterns.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/base/function_cspline.h>

#include <deal.II/base/parameter_handler.h>

using namespace dealii;

template<int dim>
class PostMagneticFluxDensity : public DataPostprocessorVector<dim> {

public:
    PostMagneticFluxDensity(): DataPostprocessorVector<dim>("B [T]", update_gradients | update_quadrature_points){}

    void evaluate_scalar_field(const DataPostprocessorInputs::Scalar<dim> &input_data,
                                       std::vector<Vector<double>> &computed_quantities) const override;

};

#endif //SOLVER_POSTMAGNETICFLUXDENSITY_H
