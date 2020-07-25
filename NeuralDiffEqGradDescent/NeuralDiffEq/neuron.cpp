#include "neuron.hpp"

Neuron::Neuron() {}

Neuron::Neuron(int inputs_num, double(*act_func_ptr)(double), double bias) {
	act_func_ptr_ = act_func_ptr;
	bias_ = bias;
	inputs_num_ = inputs_num;

	weights_arr_ptr_ = new double[inputs_num_];
	for (int i = 0; i < inputs_num_; i++)
		weights_arr_ptr_[i] = init_weight_;

	is_initialized_ = true;
}

Neuron::~Neuron() {
	delete[] weights_arr_ptr_;
}

void Neuron::init(int inputs_num, double(*act_func_ptr)(double), double bias) {
	act_func_ptr_ = act_func_ptr;
	bias_ = bias;
	inputs_num_ = inputs_num;

	weights_arr_ptr_ = new double[inputs_num_];
	for (int i = 0; i < inputs_num_; i++)
		weights_arr_ptr_[i] = init_weight_;

	is_initialized_ = true;
}

double Neuron::calc_output(const std::vector<double> &inputs_vec) {
	try {
		if (!is_initialized_)
			throw 1;

		double sum = 0;
		for (int i = 0; i < inputs_num_; i++)
			sum += weights_arr_ptr_[i] * inputs_vec[i];
		sum += bias_;

		return act_func_ptr_(sum);
	}
	catch(int e) {
		switch(e) {
			case 1:
			    std::cout << "ERROR: NEURON WAS NOT INITIALIZED" << std::endl;
			    exit(-1);
			default:
				std::cout << "ERROR: CALCULATING NEURON OUTPUT" << std::endl;
			    exit(-1);
		}
		
	}	
}

int Neuron::get_inputs_num() {
	return inputs_num_;
}

double Neuron::get_weight(int j) {
	try {
		if (!is_initialized_)
			throw 1;

		if (j > inputs_num_ || j<=0)
			throw 2;

		return weights_arr_ptr_[j-1];

	}
	catch (int e) {
		switch (e) {
		case 1:
			std::cout << "ERROR: NEURON WAS NOT INITIALIZED" << std::endl;
			exit(-1);
		case 2:
			std::cout << "ERROR: WRONG NEURON OUTPUT NUMBER" << std::endl;
			exit(-1);
		default:
			std::cout << "ERROR: CALCULATING NEURON OUTPUT" << std::endl;
			exit(-1);
		}

	}
}

void Neuron::set_weight(int j, double weight) {
	try {
		if (!is_initialized_)
			throw 1;

		weights_arr_ptr_[j-1] = weight;

	}
	catch (int e) {
		switch (e) {
		case 1:
			std::cout << "ERROR: NEURON WAS NOT INITIALIZED" << std::endl;
			exit(-1);
		default:
			std::cout << "ERROR: CALCULATING NEURON OUTPUT" << std::endl;
			exit(-1);
		}

	}

}

double Neuron::get_bias() {
	return bias_;
}

void Neuron::set_bias(double bias) {
	bias_ = bias;
}