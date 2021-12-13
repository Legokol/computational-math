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

double GaussQuadrature(double a, double b) {
    std::vector<double> w = {0.2025782419255613,
                             0.1984314853271116,
                             0.1984314853271116,
                             0.1861610000155622,
                             0.1861610000155622,
                             0.1662692058169939,
                             0.1662692058169939,
                             0.1395706779261543,
                             0.1395706779261543,
                             0.1071592204671719,
                             0.1071592204671719,
                             0.0703660474881081,
                             0.0703660474881081,
                             0.0307532419961173,
                             0.0307532419961173
    };
    std::vector<double> x = {
            0.0000000000000000,
            -0.2011940939974345,
            0.2011940939974345,
            -0.3941513470775634,
            0.3941513470775634,
            -0.5709721726085388,
            0.5709721726085388,
            -0.7244177313601701,
            0.7244177313601701,
            -0.8482065834104272,
            0.8482065834104272,
            -0.9372733924007060,
            0.9372733924007060,
            -0.9879925180204854,
            0.9879925180204854
    };

    double I = 0;
    for (int i = 0; i < x.size(); ++i) {
        I += w[i] * f((a + b) / 2 + (b - a) / 2 * x[i]);
    }
    I *= (b - a) / 2;
    return I;
}

int main() {
    double a = 0;
    double b = 3;
    std::cout << "Rectangle rule: " << rectangleRule(30000, a, b) << std::endl;
    std::cout << "Trapezoidal rule: " << trapezoidalRule(30000, a, b) << std::endl;

    int n = 15;
    double I = 0;
    double h = (b - a) / n;
    for (int i = 0; i < n; ++i) {
        I += GaussQuadrature(a + h * i, a + h * (i + 1));
    }
    std::cout << "Gaussian quadrature: " << I << std::endl;
    return 0;
}