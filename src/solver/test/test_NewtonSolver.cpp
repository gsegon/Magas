//
// Created by gordan on 3/12/23.
//
#include <gtest/gtest.h>
#include <tuple>
#include <unordered_map>
#include <string>
#include <fstream>
#include <variant>
#include <filesystem>

#include "NewtonSolver.h"
#include "LinearNuCurve.h"
#include "AnalyticNuCurve.h"
#include "ConstFSource.h"


using namespace dealii;

TEST(TestNewtonSolver, instantiation){

    NewtonSolver<2> solver;

}

TEST(TestNewtonSolver, read_mesh){
    std::filesystem::path home = std::getenv("HOME");
    std::filesystem::path test_mesh = "../../../examples/unit_square/unit_square.msh";

    NewtonSolver<2> solver;
    solver.read_mesh(test_mesh);


}

TEST(TestNewtonSolver, setup_system){

    std::filesystem::path home = std::getenv("HOME");
    std::filesystem::path test_mesh = "../../../examples/unit_square/unit_square.msh";

    NewtonSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.setup_system(true);

}

TEST(TestNewtonSolver, set_maps){

    std::filesystem::path home = std::getenv("HOME");
    std::filesystem::path test_mesh = "../../../examples/unit_square/unit_square.msh";
    std::unordered_map<int, NuCurve*> nu_map{{6, new LinearNuCurve{1}}};
    std::unordered_map<int, std::variant<FSource*, std::pair<double, double>>> f_map{{6, new ConstFSource{1}}};
    std::unordered_map<int, double> dc_map{{5, 0}};

    NewtonSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.set_nu_map(nu_map);
    solver.set_f_map(f_map);
    solver.set_dc_map(dc_map);
    solver.setup_system(true);


}

TEST(TestNewtonSolver, assemble){

    std::filesystem::path home = std::getenv("HOME");
    std::filesystem::path test_mesh = "../../../examples/unit_square/unit_square.msh";
    std::unordered_map<int, NuCurve*> nu_map{{6, new LinearNuCurve{1}}};
    std::unordered_map<int, std::variant<FSource*, std::pair<double, double>>> f_map{{6, new ConstFSource{1}}};
    std::unordered_map<int, double> dc_map{{5, 0}};

    NewtonSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.set_nu_map(nu_map);
    solver.set_f_map(f_map);
    solver.set_dc_map(dc_map);
    solver.setup_system(true);
    solver.assemble_system();

}

TEST(TestNewtonSolver, solve_initial){

    std::filesystem::path home = std::getenv("HOME");
    std::filesystem::path test_mesh = "../../../examples/unit_square/unit_square.msh";
    std::unordered_map<int, NuCurve*> nu_map{{6, new LinearNuCurve{1}}};
    std::unordered_map<int, std::variant<FSource*, std::pair<double, double>>> f_map{{6, new ConstFSource{1}}};
    std::unordered_map<int, double> dc_map{{5, 0}};

    NewtonSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.set_nu_map(nu_map);
    solver.set_f_map(f_map);
    solver.set_dc_map(dc_map);
    solver.setup_system(true);
    solver.assemble_system();
    solver.solve(0.1);

}

TEST(TestNewtonSolver, solve){

    std::filesystem::path home = std::getenv("HOME");
    std::filesystem::path test_mesh = "../../../examples/unit_square/unit_square.msh";
    std::unordered_map<int, NuCurve*> nu_map{{6, new LinearNuCurve{1}}};
    std::unordered_map<int, std::variant<FSource*, std::pair<double, double>>> f_map{{6, new ConstFSource{1}}};
    std::unordered_map<int, double> dc_map{{5, 0}};

    NewtonSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.set_nu_map(nu_map);
    solver.set_f_map(f_map);
    solver.set_dc_map(dc_map);
    solver.setup_system(true);

    int i = 0;
    double alpha = 0.0;
    while(i++ < 50){
//        solver.setup_system(false);
        solver.assemble_system();
        if (i < 5){
            alpha = 0.1;
        }
        else if (i < 15){
            alpha = 0.3;
        }
        else{
            alpha = 0.5;
        }

        solver.solve(alpha);
        std::cout << "\tResidual(" << i<< "): " << solver.compute_residual(alpha) << std::endl;
    }


}

TEST(TestNewtonSolver, EI_core){

    constexpr double mu_0 = 1.2566370621219e-6;
    constexpr double nu_0 = 1/mu_0;
    constexpr double nu_core = 1/(2500*mu_0);
    double J1 = 1*66/8.0645e-05;
    double J2 = -1*66/8.0645e-05;

    std::filesystem::path home = std::getenv("HOME");
    std::filesystem::path test_mesh = "../../../examples/EI_core/EI_core.msh";
    std::unordered_map<int, NuCurve*> nu_map{{200, new AnalyticNuCurve{}},       // Core1
                                             {201, new AnalyticNuCurve{}},       // Core2
                                             {202, new LinearNuCurve{nu_0}},       // Copper
                                             {203, new LinearNuCurve{nu_0}},       // Copper
                                             {204, new LinearNuCurve{nu_0}},       // Air
                                             {205, new LinearNuCurve{nu_0}},       // Air
    };

    std::unordered_map<int, std::variant<FSource*, std::pair<double, double>>> f_map{ {200,  new ConstFSource{0}},        // Core1
                                                                                      {201,  new ConstFSource{0}},        // Core2
                                                                                      {202,  new ConstFSource{J1}},       // Copper
                                                                                      {203, new ConstFSource{J2}},       // Copper
                                                                                      {204, new ConstFSource{0}},        // Air
                                                                                      {205, new ConstFSource{0}},        // Air
    };

    std::unordered_map<int, double> dc_map{{44, 0}} ;      // Air

    NewtonSolver<2> solver;
    solver.read_mesh(test_mesh);
    solver.set_nu_map(nu_map);
    solver.set_f_map(f_map);
    solver.set_dc_map(dc_map);
    solver.setup_system(true);

    int i = 0;
    double alpha = 0.0;
    while(i++ < 50){
//        solver.setup_system(false);
        if (i < 10)
            alpha = 0.1;
        else if (i < 15)
            alpha = 0.3;
        else
            alpha = 0.5;

        solver.assemble_system();
        solver.solve(alpha);
        double res = solver.compute_residual(alpha);
        std::cout << "\tResidual(" << i<< "): " << res << std::endl;
        if (res < 1e-6){
            std::cout << "Converged!";
            break;
        }
    }
}