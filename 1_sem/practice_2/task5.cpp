#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include "Runge_Kutta.h"

using Eigen::VectorXd;

VectorXd f(double t, const VectorXd &x) {
    double mu = 0.012277471;
    double eta = 1 - mu;
    double A = sqrt(pow(pow(x(0) + mu, 2) + x(2) * x(2), 3));
    double B = sqrt(pow(pow(x(0) - eta, 2) + x(2) * x(2), 3));

    VectorXd res(4);
    res(0) = x(1);
    res(1) = x(0) + 2 * x(3) - eta * (x(0) + mu) / A - mu * (x(0) - eta) / B;
    res(2) = x(3);
    res(3) = x(2) - 2 * x(1) - eta * x(2) / A - mu * x(2) / B;
    return res;
}

double norm(const VectorXd &x) {
    return x.cwiseAbs().maxCoeff();
}

int main() {
    VectorXd x0(4);
    x0(0) = 0.994;
    x0(1) = 0;
    x0(2) = 0;
    x0(3) = -2.00158510637908252240537862224;
    double T = 17.0652165601579625588917206249;
    double a = 0;
    double b = 5 * T;
    double h = 1e-4;
    double epsilon = 1e-3;

    Runge_Kutta<VectorXd> solver;
    auto solution = solver.DormandPrince45(a, b, x0, f, norm, epsilon, h);

    std::ofstream sol;
    sol.open("solution.txt");
    int step = solution.size() / 10000;
    for (int i = 0; i < solution.size(); i += step) {
        sol << solution[i].x;
        for (int j = 0; j < 4; ++j) {
            sol << ',' << solution[i].y[j];
        }
        sol << '\n';


    }
    sol.close();
    return 0;
}