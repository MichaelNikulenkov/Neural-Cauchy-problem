#include "euler.hpp"

void euler(double u0, double(*f)(double u, double x), double step, std::vector<double> &t_vec, std::vector<double> &result)
{
	result.resize(t_vec.size());
	result[0] = u0;

	for (int i = 1; i < t_vec.size(); i++) {
		result[i] = result[i - 1] + step * f(result[i - 1], t_vec[i - 1]);
	}

}