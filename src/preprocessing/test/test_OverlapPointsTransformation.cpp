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

    cout << endl << endl;

    vector<Vector2f> tocke1{{-1, -1},
                            {1, -1},
                            {1, 1},
                            {-1, 1}};


    ofstream tocke1_file;
    tocke1_file.open("/home/gordan/Programs/solver/scripts/tocke1.csv");
    tocke1_file << "x" << "," << "y" << endl;
    for (auto tocka : tocke1){
        tocke1_file << tocka.x() << "," << tocka.y() << endl;
    }
    tocke1_file.close();


    Rotation2D<float> rot2(M_PI/2/2);
    auto t = Translation<float, 2>(5, 0);

    vector<Vector2f> tocke2(tocke1.size());
    ofstream tocke2_file;
    tocke2_file.open("/home/gordan/Programs/solver/scripts/tocke2.csv");
    tocke2_file << "x" << "," << "y" << endl;
    for (int i=0; i < tocke1.size(); i++){
        tocke2[i] = t*rot2*tocke1[i];
        tocke2_file << tocke2[i].x() << "," << tocke2[i].y() << endl;
    }
    tocke2_file.close();


    // ----------------------------------
    Vector2f centeroid_1{0, 0};
    Vector2f centeroid_2{0, 0};

    for (auto tocka : tocke1){
        centeroid_1 += tocka;
    }
    centeroid_1 /= tocke1.size();

    for (auto tocka : tocke2){
        centeroid_2 += tocka;
    }
    centeroid_2 /= tocke2.size();

    cout << "centeroid 1: " << centeroid_1 << endl;
    cout << "centeroid 2: " << centeroid_2 << endl;

    MatrixXf A, B;
    MatrixXf Ca, Cb;
    A.resize(2, tocke1.size());
    B.resize(2, tocke2.size());
    Ca.resize(2, tocke1.size());
    Cb.resize(2, tocke2.size());

    for (int i=0; i < tocke1.size(); i++){
        A.col(i) = tocke1[i];
        B.col(i) = tocke2[i];
        Ca.col(i) = centeroid_1;
        Cb.col(i) = centeroid_2;
    }

    auto Ap = A-Ca;
    auto Bp = B-Cb;
    auto H = Ap*Bp.transpose();
    JacobiSVD<Matrix2f, ComputeFullV | ComputeFullU> svd(H);
    auto R = svd.matrixV()*svd.matrixU().transpose();
    auto t_new = Cb-R*Ca;


    cout << "R: " << endl << R << endl;
    cout << "t_new: " << endl << t_new << endl;
    cout << endl;

    cout << "A: " << endl << A << endl;
    cout << "B: " << endl << B << endl << endl;
    cout << "transform(A) = B = R*A+t_new: " << endl << R*A+t_new << endl;

}

