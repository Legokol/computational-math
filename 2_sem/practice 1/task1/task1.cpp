//
// Created by legokol on 28.02.2022.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

struct Point {
    // Структура, описывающая точку сетки
    double t;
    double x;
};

struct Solution {
    // Структура, описывающая сеточное решение
    Point point;
    double f;
};

class LeftAngle {
public:

    std::vector<std::vector<Solution>> calc(double T, double L, double tau, double h) {
        int n_t = T / tau + 1;
        int n_x = L / h + 1;
        std::vector<std::vector<Solution>> solution(n_t);

        solution[0].resize(n_x);
        for (int j = 0; j < n_x; ++j) {
            // Начальные условия
            solution[0][j] = {{0, j * h}, sin(4 * M_PI * (j * h) / L)};
        }

        for (int i = 0; i < n_t - 1; ++i) {
            solution[i + 1].resize((n_x));
            for (int j = 1; j < n_x; ++j) {
                solution[i + 1][j] = {{(i + 1) * tau, j * h},
                                      solution[i][j].f - tau / h * (solution[i][j].f - solution[i][j - 1].f)};
            }
            solution[i + 1][0] = {{(i + 1) * tau, 0}, solution[i + 1].back().f};
        }
        return solution;
    }
};

class LaxWendroff {
public:
    std::vector<std::vector<Solution>> calc(double T, double L, double tau, double h) {
        int n_t = T / tau + 1;
        int n_x = L / h + 1;
        std::vector<std::vector<Solution>> solution(n_t);

        solution[0].resize(n_x);
        for (int j = 0; j < n_x; ++j) {
            // Начальные условия
            solution[0][j] = Solution{{0, j * h}, sin(4 * M_PI * (j * h) / L)};
        }

        for (int i = 0; i < n_t - 1; ++i) {
            solution[i + 1].resize((n_x));
            for (int j = 1; j < n_x - 1; ++j) {
                solution[i + 1][j] = {{(i + 1) * tau, j * h},
                                      solution[i][j].f - tau / (2 * h) * (solution[i][j + 1].f - solution[i][j - 1].f -
                                                                          tau / h *
                                                                          (solution[i][j + 1].f - 2 * solution[i][j].f +
                                                                           solution[i][j - 1].f))};
            }
            solution[i + 1].back() = {(i + 1) * tau, (n_x - 1) * h,
                                      solution[i].back().f -
                                      tau / (2 * h) * (solution[i][0].f - solution[i][n_x - 2].f -
                                                       tau / h * (solution[i][0].f - 2 * solution[i][n_x - 1].f +
                                                                  solution[i][n_x - 2].f))};
            solution[i + 1][0] = {{(i + 1) * tau, 0}, solution[i + 1].back().f};
        }
        return solution;
    }
};

int main() {
    double T = 18; // Временной интервал
    double L = 20; // Пространственный интервал
    double h = 0.5; // Шаг по пространству
    // Числа Куранта
    double cfl1 = 0.6;
    double cfl2 = 1;
    double cfl3 = 1.01;

    LeftAngle leftAngleSolver;
    auto leftAngleSol1 = leftAngleSolver.calc(T, L, cfl1 * h, h);
    auto leftAngleSol2 = leftAngleSolver.calc(T, L, cfl2 * h, h);
    auto leftAngleSol3 = leftAngleSolver.calc(T, L, cfl3 * h, h);

    LaxWendroff laxWendroff;

    auto LaxWendroffSol1 = laxWendroff.calc(T, L, cfl1 * h, h);
    auto LaxWendroffSol2 = laxWendroff.calc(T, L, cfl2 * h, h);
    auto LaxWendroffSol3 = laxWendroff.calc(T, L, cfl3 * h, h);

    /*std::ofstream writer;
    writer.open("LaxWendroff1.01.txt");
    for (int i = 0; i < LaxWendroffSol3.size(); ++i) {
        for (int j = 0; j < LaxWendroffSol3[i].size() - 1; ++j) {
            writer << LaxWendroffSol3[i][j].f << ',';
        }
        writer << LaxWendroffSol3[i].back().f << std::endl;
    }
    writer.close();*/

    return 0;
}