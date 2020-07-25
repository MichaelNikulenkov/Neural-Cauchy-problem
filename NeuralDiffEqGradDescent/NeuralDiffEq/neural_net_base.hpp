#pragma once
#include "neuron.hpp"
#include "act_func.hpp"

class NeuralNetBase {
public:
	~NeuralNetBase();
protected:
	double(*act_func_ptr_)(double) = sigmoid;
	double(*f_act_func_ptr_)(double) = linear;
};