#pragma once
#include <math.h>


const double u0_1 = 1;
double f1(double u, double x); 
double f_der1(double u, double x);
double f_an1(double x);

const double u0_2 = 2;
double f2(double u, double x);
double f_der2(double u, double x);
double f_an2(double x);

const double u0_3 = 0;
double f3(double u, double x);
double f_der3(double u, double x);
double f_an3(double x);

const double u0_4 = 1;
double f4(double u, double x);
double f_der4(double u, double x);
double f_an4(double x);