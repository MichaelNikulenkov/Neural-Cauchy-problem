#pragma once
#include <vector>
#include "neuron.hpp"
#include "act_func.hpp"
#include "utils.hpp"

class NeuralNetFirstOrder {
public:
	NeuralNetFirstOrder(int layer_length, int t);
	~NeuralNetFirstOrder();
	void train(double u0, double(*func_ptr_)(double u, double x), double(*func_u_der_ptr_)(double u, double x), double eps, const std::vector<double> &grid_vec, int& it_count, double& loss);
	double calc_output(double x);
	double get_output_bias();
protected:
	Neuron** net_arr_ptr_ = nullptr;
	Neuron* output_neuron_ptr_ = nullptr;
	double(*act_func_ptr_)(double) = sigmoid;
	double(*act_func_first_der_ptr_)(double) = sigmoid_first_der;
	double(*act_func_sec_der_ptr_)(double) = sigmoid_second_der;

	double(*f_act_func_ptr_)(double) = linear;

	int inputs_num_ = 1;
	int layer_length_ = 1;

	double L_threshhold = 0.0000001;

	double diff_tau_ = 0.0001;
	double init_rate_ = 0.0001;
	double bias_ = 0;

	double fragment_rate_ = 0.5;
	double inv_fragment_rate_ = 2;
};