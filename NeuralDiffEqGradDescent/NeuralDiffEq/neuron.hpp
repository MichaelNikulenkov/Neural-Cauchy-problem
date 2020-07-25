#pragma once
#include <iostream>
#include <cstdio>
#include <vector>

class Neuron {
public:
	Neuron();
	Neuron(int inputs_num, double(*act_func_ptr)(double), double bias = 0);
	~Neuron();
	void init(int inputs_num, double(*act_func_ptr)(double), double bias = 0);
	double calc_output(const std::vector<double> &inputs_vec);
	int get_inputs_num();
	double get_weight(int j);
	double get_bias();
	void set_weight(int j, double weight);
	void set_bias(double bias);
private:
	double(*act_func_ptr_)(double) = nullptr;
	double *weights_arr_ptr_ = nullptr;
	double bias_ = 0;
	double init_weight_ = 1;
	
	int inputs_num_ = 0;
	bool is_initialized_ = false;
};