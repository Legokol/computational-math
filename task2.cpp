#include <iostream>
#include <cmath>

using namespace std;

const double fm = 1 / sqrt(2 * exp(1));

double f1(double x) {
    return (fm / 2 * exp(x * x));
}

double f2(double x) {
    return sqrt(log(2 * x / fm));
}

int main() {
    double epsilon = 1e-3;
    double q1 = 0.5;
    double q2 = 1 / sqrt(1 + log(4));

    double x1 = 0.5;
    double x2 = 1;

    // поиск первого корня
    double x = f1(x1);
    while (fabs(x - x1) > (1 - q1) * epsilon / 2) {
        x1 = x;
        x = f1(x1);
    }
    cout << "x1 = " << x1 << endl;

    // поиск второго корня
    x = f2(x2);
    while (fabs(x - x2) > (1 - q2) * epsilon / 2) {
        x2 = x;
        x = f2(x2);
    }
    cout << "x2 = " << x2 << endl;

    cout <<"x2 - x1 = " << x2 - x1 << endl;

    return 0;
}