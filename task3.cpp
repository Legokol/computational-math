#include <iostream>
#include <fstream>
#include "interpolation/Newton.h"
#include "interpolation/CubicSlpine.h"
#include <cmath>

using namespace std;

int main() {
    vector<double> year = {1910, 1920, 1930, 1940, 1950, 1960, 1970, 1980, 1990, 2000};
    vector<double> population = {
            92228496,
            106021537,
            123202624,
            132164569,
            151325798,
            179323175,
            203211926,
            226545805,
            248709873,
            281421906
    };

    double x = 2010;
    double y = 308745538;

    // интерполяционный многочлен в форме Ньютона
    Newton interpolator(year, population);
    double y_Newton = interpolator.interpolate(2010);

    cout << "Population from Newton polynomial: " << y_Newton << endl;
    cout << "Absolute error: " << fabs(y_Newton - y) << " , relative error: " << fabs(y_Newton - y) / y << endl;

    // кубический сплайн
    CubicSpline spline(year, population);
    double y_spline = spline.interpolate(x);

    cout << "Population from spline: " << y_spline << endl;
    cout << "Absolute error: " << fabs(y_spline - y) << " , relative error: " << fabs(y_spline - y) / y << endl;

}