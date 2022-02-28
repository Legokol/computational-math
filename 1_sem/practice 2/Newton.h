#ifndef NEWTON_METHOD_NEWTON_H
#define NEWTON_METHOD_NEWTON_H

#include <Eigen/Dense>
#include <iostream>

using Eigen::Vector2d;
using Eigen::MatrixXd;

class Newton {
public:
    Newton() {}

    double solve(double u0, const std::function<double(double)> &f, const std::function<double(double)> &df, double epsilon) {
        double u = u0 - f(u0) / df(u0);
        while (fabs(u - u0) > epsilon) {
            u0 = u;
            u = u0 - f(u0) / df(u0);
        }
        return u;
    }

    double solve(double u0, double u1, const std::function<double(double)> &f, double epsilon) {
        double u = u1 - f(u1) * (u1 - u0) / (f(u1) - f(u0));
        while (fabs(u - u1) > epsilon) {
            u0 = u1;
            u1 = u;
            u = u1 - f(u1) * (u1 - u0) / (f(u1) - f(u0));
        }
        return u;
    }

    Vector2d solve(Vector2d &u0, const std::function<Vector2d(const Vector2d &)>& f, const std::function<MatrixXd(const Vector2d &)> &J, double epsilon, bool invert = true) {
        Vector2d u;
        MatrixXd J0 = J(u0).inverse();
        u = u0 - J0 * f(u0);
        while ((u - u0).cwiseAbs().maxCoeff() > epsilon) {
            u0 = u;
            if (invert)
                u = u0 - J(u0).inverse() * f(u0);
            else
                u = u0 - J0 * f(u0);
        }
        return u;
    }
};


#endif //NEWTON_METHOD_NEWTON_H
