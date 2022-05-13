//
// Created by legokol on 14.05.22.
//

#include <fstream>
#include <iostream>
#include <vector>

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
    double p;
};

// Функция для расчёта плотности по известному давлению
double calcDensity(double pressure) {
    double rho0 = 1e3;
    double p0 = 120e5;
    double cf = 10e-4 * 1e-5;

    return rho0 * (1 + cf * (pressure - p0));
}

// Функция для аппроксимации плотности между узлами
double approximateDensity(const Node &left, const Node &right) {
    if (left.p >= right.p) {
        return calcDensity(left.p);
    } else {
        return calcDensity(right.p);
    }
}

int main() {
    double L = 500; // Длина пласта
    double k = 1e-14; // Проницаемость пласта
    double mu = 1e-3; // Вязкость жидкости
    double phi = 0.2; // Пористость пласта
    double cf = 10e-4 * 1e-5; // Сжимаемость жидкости

    double p0 = 100e5; // Начальное условие
    double pInj = 150e5; // Левое граничное условие
    double pProd = 50e5; // Правое граничное условие

    double endTime = 10 * 24 * 3600; // Время остановки в секундах

    int M = 101; // Число узлов

    double h = L / (M - 1); // Шаг по пространству
    double tau = 3600; // Шаг по времени

    double t = 0;
    std::vector<std::vector<Node>> solution;
    solution.push_back(std::vector<Node>(M));

    for (int i = 0; i < M; ++i) {
        solution[0][i] = {t, i * h, p0};
    }
    return 0;
}