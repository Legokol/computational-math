#include "Runge_Kutta.h"
#include "Newton.h"
#include <iostream>
#include <fstream>
#include <Eigen/Dense>

using Eigen::Vector2d;

Vector2d f(double t, const Vector2d &x) {
    return Vector2d(x(1), 1000 * (1 - x(1) * x(1)) * x(1) - x(0));
}

Eigen::MatrixXd J(const Vector2d &x) {
    Eigen::MatrixXd res(2, 2);
    res(0, 0) = 0;
    res(0, 1) = 1;
    res(1, 0) = -1;
    res(1, 1) = 1000 * (1 - 3 * x(1) * x(1));
    return res;
}

double norm(const Vector2d &x) {
    return x.cwiseAbs().maxCoeff();
}

/*std::vector<point<T>> GaussLegendre4(double a, double b, Vector2d &y0, Vector2d (*f)(double, Vector2d), double h, double epsilon) {
    std::array<double, 2> c = {.5 - sqrt(3) / 6, .5 + sqrt(3) / 6};
    std::array<double, 2> B = {.5, .5};
    std::array<std::array<double, 2>, 2> A = {.25, .25 - sqrt(3) / 6,
                                              .25 + sqrt(3) / 6, .25};

    std::vector<point<Vector2d>> solution((b - a) / h + 1);
    for (int i = 0; i < solution.size(); ++i) {
        solution[i].x = a + i * h;
    }
    solution[0].y = y0;
}*/

std::vector<point<Vector2d>>
backwardEuler(double a, double b, Vector2d &y0, const std::function<Vector2d(double, const Vector2d &)> &f, double h, double epsilon) {
    std::vector<point<Vector2d>> solution((b - a) / h + 1);
    for (int i = 0; i < solution.size(); ++i) {
        solution[i].x = a + i * h;
    }
    solution[0].y = y0;
    for (int i = 1; i < solution.size(); ++i) {
        Newton solver;
        Vector2d x = solution[i - 1].y + f(a, y0) * h / 2;
        auto y = solver.solve(x,
                              [&, solution, i, h](const Vector2d &x) { return x - f(solution[i].x, x) * h - solution[i - 1].y; },
                              [&h](const Vector2d &x) {
                                  MatrixXd res = MatrixXd::Identity(2, 2);
                                  res = res - h * J(x);
                                  //std::cout << res << std::endl;
                                  return res;
                              }, epsilon);
        solution[i].y = y;
    }
    return solution;
}

int main() {
    Vector2d x0(0, 1e-3);
    double a = 0;
    double b = 1000;
    double h = 1;
    double epsilon = 1e-3;
    Runge_Kutta<Vector2d> solver;

    /*auto solution1 = backwardEuler(a, b, x0, f, h*1e-1, epsilon);

    std::ofstream sol1;
    sol1.open("solution1.txt");
    for (int i = 0; i < solution1.size(); ++i) {
        sol1 << solution1[i].x << ',' << solution1[i].y(0) << '\n';
    }
    sol1.close();*/

    auto solution2 = solver.DormandPrince45(a, b, x0, f, norm, epsilon, h);


    std::ofstream sol2;
    sol2.open("solution2.txt");
    for (int i = 0; i < solution2.size(); ++i) {
        sol2 << solution2[i].x << ',' << solution2[i].y(0) << '\n';
    }
    sol2.close();
    return 0;
}