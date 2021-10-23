#include <iostream>
#include <vector>
#include <cmath>

double f(double x) {
    return (sin(100 * x) * cos(2 * x) * exp(-x * x));
}

double rectangleRule(int n, double a, double b) {
    double h = (b - a) / n;
    double I = 0;
    for (int i = 1; i <= n; ++i) {
        I += h * f(a + h * (2 * i - 1) / 2);
    }
    return I;
}

double trapezoidalRule(int n, double a, double b) {
    double h = (b - a) / n;
    double I = 0;
    for (int i = 1; i <= n; ++i) {
        I += h * (f(a + h * (i - 1)) + f(a + i * h)) / 2;
    }
    return I;
}

int main() {
    std::cout << "Rectangle rule: " << rectangleRule(100000, 0, 3) << std::endl;
    std::cout << "Trapezoidal rule: " << trapezoidalRule(40000, 0, 3) << std::endl;
    return 0;
}