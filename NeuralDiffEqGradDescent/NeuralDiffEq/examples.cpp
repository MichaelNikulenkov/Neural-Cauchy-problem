#include "examples.hpp"

double f1(double u, double x) {
	return u * u * x;
}

double f_der1(double u, double x) {
	return 2 * u * x;
}

double f_an1(double x) {
	return  2 / (2 - x * x);
}

double f2(double u, double x) {
	return -2 * u;
}


double f_der2(double u, double x) {
	return -2;
}

double f_an2(double x) {
	return 2 * exp(-2 * x);
}

double f3(double u, double x) {
	return -u / 5 + exp(-x/5) * cos(x);
}

double f_der3(double u, double x) {
	return -1 / 5;
}

double f_an3(double x) {
	return exp(-x/5) * sin(x);
}

double f4(double u, double x) {
	return x * x*x + 2 * x + x * x*(1 + 3 * x*x) / (1 + x + x * x*x) - u * (x + (1 + 3 * x * x) / (1 + x + x * x * x));
}


double f_der4(double u, double x) {
	return (1 + 3 * x * x) / (1 + x + x * x * x);
}

double f_an4(double x) {
	return x*x + exp(-x*x * 0.5) / (1 + x + x*x*x);
}