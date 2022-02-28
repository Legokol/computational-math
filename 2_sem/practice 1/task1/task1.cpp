//
// Created by legokol on 28.02.2022.
//

#include <iostream>
#include <vector>
#include <cmath>

struct Point {
    // Структура, описывающая точку сетки
    double x;
    double t;
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
        grid.resize(L / h + 1);
        std::vector<std::vector<Solution>> solution(grid.size());

        for (int i = 0; i < grid.size(); ++i) {
            // Заполнение массива сетки и начальных условий
            grid[i].resize(T / tau + 1);
            for (int j = 0; j < grid[i].size(); ++j) {
                grid[i][j].x = i * h;
                grid[i][j].t = j * tau;
            }
            solution[i].resize(grid[i].size());
            solution[i][0] = Solution{{grid[i][0].x, 0}, sin(4 * M_PI * grid[i][0].x / L)};
        }
        return solution;
    }

    //Point step() {}
};

int main() {
    return 0;
}