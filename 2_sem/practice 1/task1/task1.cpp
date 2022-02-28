//
// Created by legokol on 28.02.2022.
//

#include <iostream>
#include <vector>
#include <cmath>

struct Point {
    // Структура, описывающая точку сетки
    double t;
    double x;
};

struct Solution {
    // Структура, описывающая сеточное решение
    Point p;
    double f;
};

class LeftAngle {
private:
    // Расчётная сетка
    std::vector<std::vector<Point>> grid;

public:

    std::vector<std::vector<Solution>> calc(double T, double L, double tau, double h) {
        grid.resize(T / tau + 1);
        std::vector<std::vector<Solution>> solution(grid.size());

        grid[0].resize(L / h + 1);
        solution[0].resize(grid[0].size());
        for (int j = 0; j < grid[0].size(); ++j) {
            // Задание начальных условий
            grid[0][j].t = 0;
            grid[0][j].x = h * j;
            solution[0][j] = Solution{grid[0][j], sin(4 * M_PI * grid[0][j].x / L)};
        }

        for (int i = 0; i < grid.size() - 1; ++i) {
            grid[i + 1].resize(grid[i].size());
            solution[i + 1].resize(solution[i].size());
            for (int j = 1; j < grid[i].size(); ++j) {
                solution[i + 1][j] = Solution{grid[i + 1][j],
                                              solution[i][j].f + tau / h * (solution[i][j].f - solution[i][j - 1].f)};
            }
            solution[i + 1][0] = Solution{grid[i + 1][0], solution[i + 1].back().f};
        }

        return solution;
    }
};

int main() {
    double T = 18;
    double L = 20;
    double tau = 1;
    double h = 1;
    LeftAngle leftAngleSolver;
    auto res = leftAngleSolver.calc(T, L, tau, h);
    return 0;
}