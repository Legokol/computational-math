//
// Created by legokol on 21.03.2022.
//

#include <iostream>
#include <Eigen/Dense>
#include <vector>

struct Point {
    // Структура, описывающая точку сетки
    double t;
    double x;
};

struct Solution {
    // Структура, описывающая сеточное решение
    Point point;
    double p; // Давление
    Eigen::Vector3d w; // Вектор, содержащий плотность, её произведение на скорость и энергию
};

int main() {
    double gamma = 5. / 3.;
    double L = 10;
    double T = 0.02;

    double pL = 10e5;
    double pR = 1e5;
    double rhoL = 13;
    double rhoR = 1.3;
    double vL = 0;
    double vR = 0;

    double h = 0.1;
    double tau = 1e-3;

    std::vector<std::vector<Solution>> solution(1);

    solution[0].resize(L / h + 1);
    for (int j = 0; j < solution[0].size(); ++j) {
        solution[0][j].point = {0, h * j};
        solution[0][j].p = pL;
        solution[0][j].w = {rhoL, rhoL * vL, pL / (gamma - 1)};
    }

    double t = 0;
    while (t + tau < T) {
        std::vector<Solution> timeLayer(solution[0].size());
        for (int j = 1; j < timeLayer.size() - 1; ++j) {
            double u = solution.back()[j].w(1) / solution.back()[j].w(0); // Скорость газа
            double e = solution.back()[j].w(2) / solution.back()[j].w(0); // Энергия газа
            double c = sqrt(gamma * (gamma - 1) * e); // Скорость звука

            Eigen::DiagonalMatrix<double, 3> Lambda(u - c, u, u + c);

            Eigen::Matrix3d OmegaT{{-u * c, c,  gamma - 1},
                                   {-c * c, 0,  gamma - 1},
                                   {u * c,  -c, gamma - 1}};

            Eigen::Matrix3d A{{0,              1,         0},
                              {-u * u,         2 * u,     gamma - 1},
                              {-u * e * gamma, gamma * e, u}};

            double lambda = (Lambda.toDenseMatrix().cwiseAbs()).maxCoeff();
            if (tau * lambda / h > 1) {
                tau = 0.9 * h / lambda;
            }

            if (T - t < tau) {
                tau = T - t;
            }

            t += tau;

            timeLayer[j].point = {t, h * j};
            timeLayer[j].w = solution.back()[j].w -
                             tau * A * (solution.back()[j + 1].w - solution.back()[j - 1].w) / 2 / h +
                             tau * (OmegaT.inverse() * (Lambda.toDenseMatrix().cwiseAbs()) * OmegaT) *
                             (solution.back()[j + 1].w - 2 * solution.back()[j].w + solution.back()[j - 1].w) / 2 / h;
            timeLayer[j].p = (gamma - 1) * timeLayer[j].w(2);
        }
        timeLayer[0] = {{t, 0}, timeLayer[1].p, timeLayer[1].w};
        timeLayer.back() = {{t, (timeLayer.size() - 1) * h},
                            timeLayer[timeLayer.size() - 2].p,
                            timeLayer[timeLayer.size() - 2].w};
        solution.push_back(timeLayer);
    }

    return 0;
}