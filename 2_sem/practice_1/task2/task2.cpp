//
// Created by legokol on 21.03.2022.
//

#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <fstream>

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
    double tau = 1e-5;

    std::vector<std::vector<Solution>> solution(1);

    solution[0].resize(2 * L / h + 1);
    for (int j = 0; j < solution[0].size(); ++j) {
        solution[0][j].point = {0, -L + h * j};
        if (-L + h * j < 0) {
            solution[0][j].p = pL;
            solution[0][j].w = {rhoL, rhoL * vL, pL / (gamma - 1)};
        } else {
            solution[0][j].p = pR;
            solution[0][j].w = {rhoR, rhoR * vR, pR / (gamma - 1)};
        }
    }

    double t = 0;
    t += tau;
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

            timeLayer[j].point = {t, -L + h * j};
            timeLayer[j].w = solution.back()[j].w -
                             tau * A * (solution.back()[j + 1].w - solution.back()[j - 1].w) / 2 / h +
                             tau * (OmegaT.inverse() * (Lambda.toDenseMatrix().cwiseAbs()) * OmegaT) *
                             (solution.back()[j + 1].w - 2 * solution.back()[j].w + solution.back()[j - 1].w) / 2 / h;
            timeLayer[j].p = (gamma - 1) * timeLayer[j].w(2);
        }
        timeLayer[0] = {{t, -L}, timeLayer[1].p, timeLayer[1].w};
        timeLayer.back() = {{t, -L + (timeLayer.size() - 1) * h},
                            timeLayer[timeLayer.size() - 2].p,
                            timeLayer[timeLayer.size() - 2].w};
        solution.push_back(timeLayer);
        if (T - t < tau) {
            tau = T - t;
        }
        t += tau;
    }

//    std::vector<std::ofstream> writer(4);
//    writer[0].open("p.text");
//    writer[1].open("rho.txt");
//    writer[2].open("u.txt");
//    writer[3].open("e.txt");
//    for (int i = 0; i < solution.size(); ++i) {
//        writer[0] << solution[i][0].point.t << ' ';
//        writer[1] << solution[i][0].point.t << ' ';
//        writer[2] << solution[i][0].point.t << ' ';
//        writer[3] << solution[i][0].point.t << ' ';
//        for (int j = 0; j < solution[i].size(); ++j) {
//            writer[0] << solution[i][j].p << ' ';
//            writer[1] << solution[i][j].w(0) << ' ';
//            writer[2] << solution[i][j].w(1) / solution[i][j].w(0) << ' ';
//            writer[3] << solution[i][j].w(2) / solution[i][j].w(0) << ' ';
//        }
//        writer[0] << '\n';
//        writer[1] << '\n';
//        writer[2] << '\n';
//        writer[3] << '\n';
//    }
//    writer[0].close();
//    writer[1].close();
//    writer[2].close();
//    writer[3].close();

    return 0;
}