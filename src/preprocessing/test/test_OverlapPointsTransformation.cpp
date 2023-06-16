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
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/base/work_stream.h>
#include <deal.II/base/multithread_info.h>



using namespace Eigen;
using namespace std;
using namespace dealii;

template<typename T>
class OverLapPointsTransformation{

public:
    OverLapPointsTransformation(std::vector<T> dataset_a, std::vector<T> dataset_b){
        // Calculates R and t such that RA+t = B. So, for A->B.

        firsts = dataset_a;
        seconds = dataset_b;

        auto centeroid_firsts = calc_centeriod(firsts);
        auto centeroid_seconds = calc_centeriod(seconds);

        A.resize(2, firsts.size());
        B.resize(2, firsts.size());
        Ca.resize(2, seconds.size());
        Cb.resize(2, seconds.size());

        for (int i=0; i < firsts.size(); i++){
            A.col(i) = (Vector2d){firsts[i][0], firsts[i][1]};
            B.col(i) = (Vector2d){seconds[i][0], seconds[i][1]};
            Ca.col(i) = centeroid_firsts;
            Cb.col(i) = centeroid_seconds;
        }

        auto Ap = A-Ca;
        auto Bp = B-Cb;
        auto H = Ap*Bp.transpose();
        JacobiSVD<Matrix2d, ComputeFullV | ComputeFullU> svd(H);
        R = svd.matrixV()*svd.matrixU().transpose();
        t = Cb-R*Ca;

    }

    Vector2d calc_centeriod(std::vector<T> points){
        Vector2d centeroid{0,0};
        for (auto point : points){
            centeroid[0] += point[0];
            centeroid[1] += point[1];
        }
        centeroid[0] /= points.size();
        centeroid[1] /= points.size();

        return centeroid;
    }

    void apply_transform(std::vector<T>& a_transformed){
        auto A_transformed = R*A+t;
        for (int i=0; i<A_transformed.cols(); i++)
            a_transformed[i] = {A_transformed.col(i)[0],  A_transformed.col(i)[1]};
    }


private:
    std::vector<T> firsts;
    std::vector<T> seconds;

    MatrixXd A, B;
    MatrixXd Ca, Cb;

    Matrix2d R;
    Matrix2Xd t;
};


TEST(OverlapPointsTransformation, basic){


    vector<Vector2f> tocke1{{-1, -1},
                            {1, -1},
                            {1, 1},
                            {-1, 1}};

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

    std::vector<temp_def> a_trans(tocke1.size());
    OverLapPointsTransformation<temp_def> olpt{a, b};
    olpt.apply_transform(a_trans);

    for (auto val: b)
        cout << val << endl;
    for (auto val: a_trans)
        cout << val << endl;

}

