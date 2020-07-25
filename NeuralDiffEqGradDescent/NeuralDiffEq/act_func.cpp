#include "act_func.hpp"

double sigmoid(double x) {
	return (1 / (1 + exp(-x)));
}

double sigmoid_first_der(double x) {
	return exp(-x) / ((1 + exp(-x)) * (1 + exp(-x)));
}

double sigmoid_second_der(double x) {
	return ( 2 * exp(-2 * x) ) / ((1 + exp(-x)) * (1 + exp(-x)) * (1 + exp(-x))) - exp(-x) / ((1 + exp(-x)) * (1 + exp(-x)));
}

double sigmoid_wx_w_der(double x, double w) {
	return (exp(-w * x) * x) / ((1 + exp(-x*w)) * (1 + exp(-x * w)));
}

double sigmoid_der_wx_w_der(double x, double w) {
	return exp(-w * x) / ((1 + exp(-x * w)) * (1 + exp(-x * w))) 
		+ (2 * exp(-2 * w * x) * w * x) / ((1 + exp(-x * w)) * (1 + exp(-x * w))* (1 + exp(-x * w)))
		- (exp(-w * x)*w*x) / ((1 + exp(-x * w)) * (1 + exp(-x * w)));
}

double linear(double x) {
	return x;
}