//
// Created by gordan on 5/12/23.
//

#include "LinearSolver.h"

template<int dim>
void print_triangulation_info(const Triangulation<dim> &triangulation)
{

    std::cout   << "Mesh info: " << std::endl
                << " dimension: " << dim << std::endl
                << " no. of cells: " << triangulation.n_active_cells() << std::endl;
    {
        std::map<types::boundary_id, unsigned int> boundary_count;
        for (const auto &face : triangulation.active_face_iterators())
            if (face->at_boundary())
                boundary_count[face->boundary_id()]++;

        std::cout << " boundary indicators: ";
        for (const std::pair<const types::boundary_id, unsigned int> &pair : boundary_count)
        {
            std::cout << pair.first << "(" << pair.second << " times)";
        }
        std::cout << std::endl;
    }
}


int main(int argc, char* argv[]){

    constexpr double mu_0 = 1.2566370621219e-6;
    constexpr double nu_0 = 1/mu_0;
    double i_current = 1e3;
    double Jdensity = i_current / (std::pow(0.1,2) * M_PI);


    std::string test_mesh = "/home/gordan/Programs/solver/test/test_data/test_2_conductors/2_conductors_x_dense.msh";

    std::unordered_map<int, double> nu_map{{3, nu_0},       // Air
                                           {1, nu_0},       // Conductor 1
                                           {2, nu_0}};      // Conductor 2

    std::unordered_map<int, double> f_map{{3, 0},           // Air
                                          {1, Jdensity},    // Conductor 1
                                          {2, -Jdensity}};  // Conductor 2

    std::unordered_map<int, double> dc_map{{100, 0}};       // Outerbounds

    LinearSolver<2> solver;
    solver.read_mesh(test_mesh);
    print_triangulation_info(solver.get_triangulation());
    solver.setup_system();
    solver.set_nu_map(nu_map);
    solver.set_f_map(f_map);
    solver.set_dc_map(dc_map);
    solver.assemble_system();
    solver.solve();
    solver.output_results("2_conductors");


    return 0;
}