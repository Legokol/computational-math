#ifndef INTERPOLATOR_NEWTON_H
#define INTERPOLATOR_NEWTON_H

#include <vector>

using std::vector;


class Newton {
public:

    Newton(const vector<double> &x, const vector<double> &f) : x0(x) {
        if (x.size() != f.size()) {
            throw -1;
        }
        for (int i = 0; i < x.size() - 1; ++ i) {
            for (int j = i + 1; j < x.size(); ++j) {
                if (x[i] == x[j]) {
                    throw -1;
                }
            }
        }

        vector<vector<double>> dd(x0.size());
        a.resize(x0.size());
        a[0] = f[0];
        dd[0] = f;
        for (int k = 1; k < x0.size(); ++k) {
            dd[k].resize(dd[k - 1].size() - 1);
            for (int i = 0; i < dd[k].size(); ++i) {
                dd[k][i] = (dd[k - 1][i+1] - dd[k - 1][i]) / (x0[i + k] - x0[i]);
            }
            a[k] = dd[k][0];
        }
    }

    double interpolate(double x) {
        double f = a[0];
        for (int k = 1; k < a.size(); ++k) {
            double d = a[k];
            for (int i = 0; i < k; ++i) {
                d *= (x - x0[i]);
            }
            f += d;
        }
        return f;
    }

private:
    vector<double> x0;
    vector<double> a;
};


#endif //INTERPOLATOR_NEWTON_H
