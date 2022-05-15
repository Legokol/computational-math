//
// Created by legokol on 14.05.22.
//

#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>

std::vector<double>
solveTriagonalSlae(const std::vector<double> &a, const std::vector<double> &b, const std::vector<double> &c,
                   const std::vector<double> &d) {
    std::vector<double> p(d.size() - 1);
    std::vector<double> q(d.size() - 1);
    std::vector<double> x(d.size());

    p[0] = -c[0] / b[0];
    q[0] = d[0] / b[0];
    for (int i = 1; i < p.size(); ++i) {
        p[i] = -c[i] / (a[i] * p[i - 1] + b[i]);
        q[i] = (d[i] - a[i] * q[i - 1]) / (a[i] * p[i - 1] + b[i]);
    }
    x.back() = (d.back() - a.back() * q.back()) / (a.back() * p.back() + b.back());
    for (int i = x.size() - 2; i >= 0; --i) {
        x[i] = p[i] * x[i + 1] + q[i];
    }
    return x;
}

// Структура, описывающая решение
struct Node {
    double t;
    double x;
    double T;
};

// Расчёт коэффициента температуропроводности
double calcThermalDiffusivity(double temperature) {
    double c1 = 1e-2;
    double c2 = 1e-5;
    return c1 + c2 * sqrt(temperature);
}

// Аппроксимация коэффициента между узлами
double approxThermalDiffusivity(const Node &left, const Node &right) {
    return calcThermalDiffusivity((left.T + right.T) / 2);
}

int main() {
    double L = 1; // Длина
    double endTime = 100; // Время остановки

    double T0 = 273; // Начальная температура
    double TLeft = 400; // Левое граничное условие
    double TRight = 250; // Правое граничное условие

    int M = 101; // Число узлов

    double h = L / (M - 1); // Шаг по пространству
    double tau = 0.01; // Шаг по времени

    double t = 0;
    std::vector<std::vector<Node>> solution;
    solution.push_back(std::vector<Node>(M));

    for (int i = 0; i < M; ++i) {
        solution[0][i] = {t, i * h, T0};
    }

    solution[0][0].T = TLeft;
    solution[0].back().T = TRight;

    t += tau;
    while (t < endTime) {
        // Заполнение матрицы СЛАУ
        std::vector<double> a(M);
        std::vector<double> b(M);
        std::vector<double> c(M);
        std::vector<double> d(M);

        b[0] = 1;
        c[0] = 0;
        d[0] = TLeft;

        b.back() = 1;
        a.back() = 0;
        d.back() = TRight;

        for (int i = 1; i < M - 1; ++i) {
            double aLeft = approxThermalDiffusivity(solution.back()[i - 1], solution.back()[i]);
            double aRight = approxThermalDiffusivity(solution.back()[i], solution.back()[i + 1]);

            a[i] = aLeft / (h * h);
            c[i] = aRight / (h * h);
            b[i] = -a[i] - c[i] - 1 / tau;
            d[i] = -1 / tau * solution.back()[i].T;
        }

        std::vector<double> temperature = solveTriagonalSlae(a, b, c, d);
        std::vector<Node> timeLayer(temperature.size());

        for (int i = 0; i < M; ++i) {
            timeLayer[i] = {t, i * h, temperature[i]};
        }

        solution.push_back(timeLayer);
        t += tau;
    }

    std::ofstream writer;
    writer.open("temperature.txt");
    for (int i = 0; i < solution.size(); i += 50) {
        writer << solution[i][0].t << ' ';
        for (int j = 0; j < solution[i].size(); ++j) {
            writer << solution[i][j].T << ' ';
        }
        writer << '\n';
    }

    return 0;
}