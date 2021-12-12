#include <iostream>
#include <cmath>
#include <array>
#include <Eigen/Dense>


using namespace std;

array<double, 2> f(const array<double, 2> &u) {
    array<double, 2> res;
    res[0] = u[0] * u[0] + u[1] * u[1] - 1;
    res[1] = u[1] - tan(u[0]);
    return res;
}

array<double, 2> iterate(const array<double, 2> &u0, double epsilon) {
    array<double, 2> u = u0;
    array<double, 2> u_next = {atan(u[1]), sqrt(1 - u[0] * u[0])};
    while (fabs(u_next[0] - u[0]) > epsilon || fabs(u_next[1] - u[1]) > epsilon) {
        u = u_next;
        u_next = {atan(u[1]), sqrt(1 - u[0] * u[0])};
    }

    return u_next;
}

Eigen::Vector2d Newton(const Eigen::Vector2d &u) {
    Eigen::Matrix2d J;
    J(0, 0) = 2 * u(0);
    J(0, 1) = 2 * u(1);
    J(1, 0) = -1 / (1 + u(0) * u(0));
    J(1, 1) = 1;
    return (u - J.inverse() * Eigen::Vector2d(u(0) * u(0) + u(1) * u(1) - 1, u(1) - tan(u(0))));
}

int main() {
    double epsilon = 1e-6;
    array<double, 2> u0 = {0.5, 0.5};
    Eigen::Vector2d u0_(u0[0], u0[1]);
    //auto u = f(u0);
    auto u = iterate(u0, epsilon);

    cout << "(x1, y1) = " << u[0] << ", " << u[1] << endl;
    cout << "(x2, y2) = " << -u[0] << ", " << -u[1] << endl;
    auto u_ = Newton(u0_);
    while ((u_ - u0_).norm() > epsilon) {
        u0_ = u_;
        u_ = Newton(u0_);
    }
    cout << u_ << endl;

    return 0;
}