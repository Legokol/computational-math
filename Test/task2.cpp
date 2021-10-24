#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

double f(double x) {
    return (pow(x, 3) + 5 * pow(x, 2) + 10 * x + 4.5);
}

int main() {
    vector<double> x = {1.5 - 0.7*sqrt(15), 1.5, 1.5 + 0.7*sqrt(15)};
    vector<double> c = {5. / 9, 8. / 9, 5. / 9};
    double I = 0;
    for (int i = 0; i < x.size(); i++) {
        I += c[i] * f(x[i]);
    }
    cout << 3.5 * I << endl;
    return 0;
}