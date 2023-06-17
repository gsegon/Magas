//
// Created by gordan on 6/16/23.
//

#ifndef SOLVER_OVERLAPPOINTSTRANSFORMATION_H
#define SOLVER_OVERLAPPOINTSTRANSFORMATION_H

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Core>

//using namespace Eigen;

template<typename T>
class OverLapPointsTransformation{

public:
    OverLapPointsTransformation(const std::vector<T> dataset_a, const std::vector<T> dataset_b){
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
            A.col(i) = (Eigen::Vector2d){firsts[i][0], firsts[i][1]};
            B.col(i) = (Eigen::Vector2d){seconds[i][0], seconds[i][1]};
            Ca.col(i) = centeroid_firsts;
            Cb.col(i) = centeroid_seconds;
        }

        auto Ap = A-Ca;
        auto Bp = B-Cb;
        auto H = Ap*Bp.transpose();
        Eigen::JacobiSVD<Eigen::Matrix2d, Eigen::ComputeFullV | Eigen::ComputeFullU> svd(H);
        R = svd.matrixV()*svd.matrixU().transpose();
        t = Cb-R*Ca;

    }

    Eigen::Vector2d calc_centeriod(std::vector<T> points){
        Eigen::Vector2d centeroid{0,0};
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

        T point;
        for (int i=0; i<A_transformed.cols(); i++){
            point[0] = A_transformed.col(i)[0];
            point[1] = A_transformed.col(i)[1];
            a_transformed[i] = point;
        }

    }


private:
    std::vector<T> firsts;
    std::vector<T> seconds;

    Eigen::MatrixXd A, B;
    Eigen::MatrixXd Ca, Cb;

    Eigen::Matrix2d R;
    Eigen::Matrix2Xd t;
};


#endif //SOLVER_OVERLAPPOINTSTRANSFORMATION_H
