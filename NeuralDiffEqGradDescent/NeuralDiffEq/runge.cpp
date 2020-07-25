#include "runge.hpp"

float dydx(float x, float y)
{
	return((x - y) / 2);
}

void runge_kutta4(double u0, double(*f)(double u, double x), double step, std::vector<double> &t_vec, std::vector<double> &result) {

	result.resize(t_vec.size());
	result[0] = u0;
	double x0 = t_vec[0];
	double y = u0;

	for (int i = 1; i < t_vec.size(); i++) {
		double k1 = step * f(y, x0);
		double k2 = step * f(y + 0.5*k1, x0 + 0.5*step);
		double k3 = step * f(y + 0.5*k2, x0 + 0.5*step);
		double k4 = step * f(y + k3, x0 + step);

		y += (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
		x0 = t_vec[i];

		result[i] = y;
	}
}