//
// Created by gordan on 5/14/23.
//


#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/data_out.h>


using namespace dealii;

#include "PostMagneticFluxDensity.h"

template class PostMagneticFluxDensity<2>;

template<int dim>
void PostMagneticFluxDensity<dim>::evaluate_scalar_field(const DataPostprocessorInputs::Scalar<dim> &input_data,
                                                         std::vector<Vector<double>> &computed_quantities) const {

        for (unsigned int p = 0; p < input_data.solution_gradients.size(); p++){

            computed_quantities[p][0] = input_data.solution_gradients[p][1];
            computed_quantities[p][1] = -input_data.solution_gradients[p][0];

        }
}
