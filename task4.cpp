#include <iostream>
#include <cmath>
#include <array>

using namespace std;

array<double, 2> f(const array<double, 2> &u) {
    array<double, 2> res;
    res[0] = u[0] * u[0] + u[1] * u[1] - 1;
    res[1] = u[1] - tan(u[0]);
    return res;
}

array<double, 2> iterate(const array<double, 2> &u0, double epsilon) {
    array<double, 2> u = u0;
    array<double, 2> u_next = {atan(u[1]), sqrt(1 - u[0] * u[0])};
    while (fabs(u_next[0] - u[0]) > epsilon || fabs(u_next[1] - u[1]) > epsilon) {
        u = u_next;
        u_next = {atan(u[1]), sqrt(1 - u[0] * u[0])};
    }

    return u_next;
}

int main() {
    double epsilon = 1e-6;
    array<double, 2> u0 = {0.5, 0.5};
    //auto u = f(u0);
    auto u = iterate(u0, epsilon);

    cout << "(x1, y1) = " << u[0] << ", " << u[1] << endl;
    cout << "(x2, y2) = " << -u[0] << ", " << -u[1] << endl;
    return 0;
}