#pragma once
#include <math.h> 

double sigmoid(double x);
double sigmoid_first_der(double x);
double sigmoid_second_der(double x);
double sigmoid_wx_w_der(double x, double w);
double sigmoid_der_wx_w_der(double x, double w);
double linear(double x);

