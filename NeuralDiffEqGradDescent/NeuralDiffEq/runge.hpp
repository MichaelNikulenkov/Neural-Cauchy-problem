#pragma once
#include <vector>

float dydx(float x, float y);

void runge_kutta4(double u0, double(*f)(double u, double x), double step, std::vector<double> &t_vec, std::vector<double> &result);
