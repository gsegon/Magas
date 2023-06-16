//
// Created by gordan on 3/12/23.
//
#include <gtest/gtest.h>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <fstream>

using namespace Eigen;
using namespace std;

template<typename T>
class OverLapPointsTransformation{

public:
    OverLapPointsTransformation(std::vector<T> firsts, std::vector<T> seconds){
        this->firsts = firsts;
        this->seconds = seconds;
    }

private:
    std::vector<T> firsts;
    std::vector<T> seconds;
};


TEST(OverlapPointsTransformation, basic){

    std::vector<std::vector<double>> firsts{{0, 0},
                                           {0, 1},
                                           {0.3, 0.2},
                                           {0.5, 0.5},
                                           {1, 1}};

    std::vector<std::vector<double>> seconds{{0, 1},
                                            {0.3, 0.2},
                                            {0, 0},
                                            {0.5, 0.5},
                                            {1, 1}};

//    MatrixXf m = MatrixXf::Random(3,2);
//    cout << "Here is the matrix m:" << endl << m << endl;
//    JacobiSVD<MatrixXf, ComputeThinU | ComputeThinV> svd(m);
//    cout << "Its singular values are:" << endl << svd.singularValues() << endl;
//    cout << "Its left singular vectors are the columns of the thin U matrix:" << endl << svd.matrixU() << endl;
//    cout << "Its right singular vectors are the columns of the thin V matrix:" << endl << svd.matrixV() << endl;
//    Vector3f rhs(1, 0, 0);
//    cout << "Now consider this rhs vector:" << endl << rhs << endl;
//    cout << "A least-squares solution of m*x = rhs is:" << endl << svd.solve(rhs) << endl;

    cout << endl << endl;

    vector<Vector2f> tocke1{{-0.5, -0.5},
                            {0.0, -1.0},
                            {0.5, 0.8},
                            {2.5, 1.2}};

    cout << "tocke1.size(): " << tocke1.size() << endl;
    cout << "tocke1:" << endl;
    ofstream tocke1_file;
    tocke1_file.open("/home/gordan/Programs/solver/scripts/tocke1.csv");
    tocke1_file << "x" << "," << "y" << endl;
    for (auto tocka : tocke1){
        cout << tocka << endl << endl;
        tocke1_file << tocka.x() << "," << tocka.y() << endl;
    }
    tocke1_file.close();


    Rotation2D<float> rot2(M_PI/4.55);
    auto t = Translation<float, 2>(4.3, 6.4);

    cout << endl;
    cout << "tocke2:" << endl;
    cout << "tocke1.size(): " << tocke1.size() << endl;
    vector<Vector2f> tocke2(tocke1.size());
    ofstream tocke2_file;
    tocke2_file.open("/home/gordan/Programs/solver/scripts/tocke2.csv");
    tocke2_file << "x" << "," << "y" << endl;
    for (int i=0; i < tocke1.size(); i++){
        cout << "i: " << i << endl;
        tocke2[i] = rot2*t*tocke1[i];
        cout << tocke2[i] << endl << endl;
        tocke2_file << tocke2[i].x() << "," << tocke2[i].y() << endl;
    }
    tocke2_file.close();
}

