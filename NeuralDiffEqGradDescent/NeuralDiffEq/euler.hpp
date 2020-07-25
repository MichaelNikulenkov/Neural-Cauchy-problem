#pragma once
#include <vector>
void euler(double u0, double(*f)(double u, double x), double step, std::vector<double> &t_vec, std::vector<double> &result);