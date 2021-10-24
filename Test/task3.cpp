#include <iostream>
#include <vector>

using namespace std;

int main() {
    double x0 = 1.1;
    vector<double> x = {0.1, 0.55, 1.3, 1.98, 2.64};
    vector<double> f = {-2.14, -0.34, -1.7, 0.2, 1.68};

    vector<double> a(x.size());

    vector<vector<double>> dd(x.size());
    a[0] = f[0];
    dd[0] = f;
    for (int k = 1; k < x.size(); ++k) {
        dd[k].resize(dd[k - 1].size() - 1);
        for (int i = 0; i < dd[k].size(); ++i) {
            dd[k][i] = (dd[k - 1][i+1] - dd[k - 1][i]) / (x[i + k] - x[i]);
        }
        a[k] = dd[k][0];
    }

    cout << "Divided differences: ";
    for (int i = 0; i < a.size(); ++i) {
        cout << a[i] << " ";
    }
    cout << endl;

    double res = a[1] + a[2] * ((x0 - x[0]) + (x0 - x[1])) + 
    a[3] * ((x0 - x[2]) * (x0 - x[1]) + (x0 - x[2]) * (x0 - x[0]) + (x0 - x[1]) * (x0 - x[0])) +
    a[4] * ((x0 - x[3]) * (x0 - x[2]) * (x0 - x[1]) + (x0 - x[3]) * (x0 - x[2]) * (x0 - x[0]) + (x0 - x[3]) * (x0 - x[1]) * (x0 - x[0]) + (x0 - x[2]) * (x0 - x[1]) * (x0 - x[0]));

    cout << "Derivative at x0: " << res << endl;
    return 0;
}