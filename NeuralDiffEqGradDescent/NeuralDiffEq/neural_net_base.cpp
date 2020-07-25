#include "neural_net_base.hpp"

NeuralNetBase::NeuralNetBase() {
	delete output_neuron_ptr_;
}

NeuralNetBase::~NeuralNetBase() {}

//double NeuralNetBase::calc_output() {
//	try {
//		if (!is_initialized_)
//			throw 1;
//
//
//	}
//	catch (int e) {
//		switch (e) {
//		case 1:
//			std::cout << "ERROR: NEURAL NET WAS NOT INITIALIZED" << std::endl;
//			exit(-1);
//		default:
//			std::cout << "ERROR: CALCULATING NEURAL NET OUTPUT" << std::endl;
//			exit(-1);
//		}
//
//	}
//}
