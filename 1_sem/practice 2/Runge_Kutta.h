#ifndef CLASSICAL_RUNGE_KUTTA_METHOD_RUNGE_KUTTA_H
#define CLASSICAL_RUNGE_KUTTA_METHOD_RUNGE_KUTTA_H

#include <vector>
#include <array>
#include <functional>
#include <cmath>

template<typename T>
struct point {
    double x;
    T y;
};

template<typename T>
class Runge_Kutta {
private:

public:
    Runge_Kutta() {}

    std::vector<point<T>> classicalRungeKutta(double a, double b, T y0, const std::function<T(double, const T &)> &f, double h) {
        std::vector<point<T>> solution((b - a) / h + 1);
        for (int i = 0; i < solution.size(); ++i) {
            solution[i].x = a + i * h;
        }
        solution[0].y = y0;
        std::array<T, 4> k;
        for (int i = 0; i < solution.size() - 1; ++i) {
            k[0] = f(solution[i].x, solution[i].y);
            k[1] = f(solution[i].x + h / 2, solution[i].y + h * k[0] / 2);
            k[2] = f(solution[i].x + h / 2, solution[i].y + h * k[1] / 2);
            k[3] = f(solution[i].x + h, solution[i].y + h * k[3]);

            solution[i + 1].y = solution[i].y + h / 6 * (k[0] + 2 * k[1] + 2 * k[2] + k[3]);
        }
        return solution;
    }

    std::vector<point<T>> DormandPrince45(double a, double b, T y0, const std::function<T(double, const T &)> &f,
                                          const std::function<double(const T &)> &norm, double epsilon, double h0) {
        std::array<std::array<double, 6>, 6> A = {1. / 5, 0, 0, 0, 0, 0,
                                                  3. / 40, 9. / 40, 0, 0, 0, 0,
                                                  44. / 45, -56. / 15, 32. / 9, 0, 0, 0,
                                                  19372. / 6561, -25360. / 2187, 64448. / 6561, -212. / 729, 0, 0,
                                                  9017. / 3168, -355. / 33, 46732. / 5247, 49. / 176, -5103. / 18656, 0,
                                                  35. / 384, 0, 500. / 1113, 125. / 192, -2187. / 6784, 11. / 84};
        std::array<double, 6> c = {1. / 5, 3. / 10, 4. / 5, 8. / 9, 1, 1};
        std::array<double, 7> b1 = {35. / 384, 0, 500. / 1113, 125. / 192, -2187. / 6784, 11. / 84, 0};
        std::array<double, 7> b2 = {5179. / 57600, 0, 7571. / 16695, 393. / 640, -92097. / 339200, 187. / 2100, 1. / 40};
        std::array<T, 7> k;
        double h = h0;
        double x = a;
        std::vector<point<T>> solution;
        solution.push_back(point<T>{a, y0});

        while (x < b) {
            k[0] = f(x, solution.back().y);
            T y1 = solution.back().y + h * b1[0] * k[0];
            T y2 = solution.back().y + h * b2[0] * k[0];
            for (int i = 1; i < 7; ++i) {
                T y = solution.back().y;
                for (int j = 0; j < i; ++j) {
                    y += h * A[i - 1][j] * k[j];
                }
                k[i] = f(x + c[i - 1] * h, y);
                y1 += h * b1[i] * k[i];
                y2 += h * b2[i] * k[i];
            }
            double err = norm(y1 - y2);
            double h_opt = 0.9 * h * pow(epsilon / err, 0.2);
            if (h_opt < h / 2)
                h = h_opt;
            else {
                x += h;
                solution.push_back(point<T>{x, y1});
                h = std::min(h_opt, h0);
                if (x + h > b) {
                    h = b - x;
                }
            }
        }
        return solution;
    }
};

#endif //CLASSICAL_RUNGE_KUTTA_METHOD_RUNGE_KUTTA_H
