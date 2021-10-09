#include <iostream>
#include <cmath>

using namespace std;

double f(double x) {
    return (x * x + pow(tan(x), 2) - 1);
}

double secantMethod(double (*f)(double), double a, double b, double epsilon) {
    if (f(a) * f(b) > 0)
        throw -1;
    double c = (a * f(b) - b * f(a)) / (f(b) - f(a));
    while (b - a > epsilon) {

        if (f(a) * f(c) < 0) {
            b = c;
            c = (a * f(b) - b * f(a)) / (f(b) - f(a));
        } else if (f(b) * f(c) < 0){
            a = c;
            c = (a * f(b) - b * f(a)) / (f(b) - f(a));
        }
        else break;
    }
    return c;
}

int main() {
    double a = 0;
    double b = 0.7;
    double epsilon = 1e-6;

    double x = secantMethod(f, a, b, epsilon);
    double y = tan(x);
    cout << "(x1, y1) = " << x << ", " << y << endl;
    cout << "(x2, y2) = " << -x << ", " << -y << endl;
    return 0;
}