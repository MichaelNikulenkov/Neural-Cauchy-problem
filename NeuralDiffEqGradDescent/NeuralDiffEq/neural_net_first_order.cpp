#include "neural_net_first_order.hpp"

NeuralNetFirstOrder::NeuralNetFirstOrder(int layer_length, int t) {
	init(t, bias_);
	layer_length_ = layer_length;
	net_arr_ptr_ = new Neuron*[layer_length_];

	for (int i = 0; i < layer_length_; i++)
		net_arr_ptr_[i] = new Neuron(inputs_num_, act_func_ptr_, 0);

	output_neuron_ptr_ = new Neuron(layer_length_, f_act_func_ptr_, bias_);


}

double NeuralNetFirstOrder::get_output_bias() {
	return output_neuron_ptr_->get_bias();
}

NeuralNetFirstOrder::~NeuralNetFirstOrder() {
	for (int i = 0; i < layer_length_; i++)
		delete net_arr_ptr_[i];

	delete[] net_arr_ptr_;
	delete output_neuron_ptr_;
}

double NeuralNetFirstOrder::calc_output(double x) {
	
	double sum = 0;
	std::vector<double> z;

	for (int i = 0; i < layer_length_; i++)
		z.push_back(net_arr_ptr_[i]->get_weight(1) * x + net_arr_ptr_[i]->get_bias());

	for (int i = 0; i < layer_length_; i++)
		sum += act_func_ptr_(z[i]) * output_neuron_ptr_->get_weight(i + 1) + output_neuron_ptr_->get_bias();

	return sum;
}

void NeuralNetFirstOrder::train(double u0, double(*func_ptr_)(double u, double x), double(*func_u_der_ptr_)(double u, double x), double eps, const std::vector<double> &grid_vec, int& it_count, double& loss) {

	// w - u - v layout
	std::vector<double> param_vec(layer_length_*3);
	std::vector<double> prev_param_vec(layer_length_ * 3, 0.0);
	std::vector<double> grad_vec(layer_length_ * 3);

	double g = output_neuron_ptr_->get_bias(); //output neuron bias

	//store output values from hidden layer z = wi*x+ui
	std::vector<std::vector<double>> z(grid_vec.size());

	//store net outputs
	std::vector<double> N (grid_vec.size(), 0.0);

	//store f(x,N)
	std::vector<double> f (grid_vec.size());

	//store dfdu
	std::vector<double> dfdu(grid_vec.size());

	//store dNdx
	std::vector<double> dNdx(grid_vec.size(), 0.0);

	//store gradient matrices and net derivatives matrices
	std::vector<std::vector<double>> dEdw(grid_vec.size());
	std::vector<std::vector<double>> dEdv(grid_vec.size());
	std::vector<std::vector<double>> dEdu(grid_vec.size());

	std::vector<std::vector<double>> dNdw(grid_vec.size());
	std::vector<std::vector<double>> dNdv(grid_vec.size());
	std::vector<std::vector<double>> dNdu(grid_vec.size());

	std::vector<std::vector<double>> dNdwdx(grid_vec.size());
	std::vector<std::vector<double>> dNdvdx(grid_vec.size());
	std::vector<std::vector<double>> dNdudx(grid_vec.size());


	for (int i = 0; i < grid_vec.size(); i++) {
		z[i].resize(layer_length_);
		dEdw[i].resize(layer_length_);
		dEdv[i].resize(layer_length_);
		dEdu[i].resize(layer_length_);
		dNdw[i].resize(layer_length_);
		dNdv[i].resize(layer_length_);
		dNdu[i].resize(layer_length_);
		dNdwdx[i].resize(layer_length_);
		dNdvdx[i].resize(layer_length_);
		dNdudx[i].resize(layer_length_);
	}

	double rate = init_rate_;
	bool is_first_iter = true;
	double L = 0;
	double grad = 0;

	//store weights and biases in vectors
	for (int i = 0; i < layer_length_; i++) {
		//w
		param_vec[i] = net_arr_ptr_[i]->get_weight(1);
		//u
		param_vec[i + layer_length_] = net_arr_ptr_[i]->get_bias();
		//v
		param_vec[i + 2 * layer_length_] = output_neuron_ptr_->get_weight(i + 1);
	}


	bool stop_flag = false;
	it_count = 0;



	while (!stop_flag) {

		//clear data

		std::fill(grad_vec.begin(), grad_vec.end(), 0.0);


		//---------------------------------------------------------------
		

		if (is_first_iter) {

			//calculate z
			for (int j = 0; j < grid_vec.size(); j++)
				for (int i = 0; i < layer_length_; i++) {
					z[j][i] = param_vec[i] * grid_vec[j] + param_vec[i + layer_length_];
				}


			//calculate net outputs
			for (int j = 0; j < grid_vec.size(); j++) {
				for (int i = 0; i < layer_length_; i++) {
					N[j] += act_func_ptr_(z[j][i]) * param_vec[i + 2 * layer_length_];
				}

				N[j] += g;
				N[j] = f_act_func_ptr_(N[j]);
			}

			//calculate dNdx
			for (int j = 0; j < grid_vec.size(); j++)
				for (int i = 0; i < layer_length_; i++) {
					dNdx[j] += param_vec[i + 2 * layer_length_] * param_vec[i] * act_func_first_der_ptr_(z[j][i]);
				}

			//calculate f(x, N)
			for (int i = 0; i < grid_vec.size(); i++) {
				f[i] = func_ptr_(u0 + grid_vec[i] * N[i], grid_vec[i]);
			}
		}

		//calculate dfdu
		for (int i = 0; i < grid_vec.size(); i++) {
			dfdu[i] = func_u_der_ptr_(u0 + grid_vec[i] * N[i], grid_vec[i]);
		}

		//calculate dNdw, dNdv, dNdu
		for (int j = 0; j < grid_vec.size(); j++)
			for (int i = 0; i < layer_length_; i++) {
				dNdw[j][i] = grid_vec[j] * param_vec[i + 2 * layer_length_] * act_func_first_der_ptr_(z[j][i]);
				dNdv[j][i] = grid_vec[j] * act_func_ptr_(z[j][i]);
				dNdu[j][i] = param_vec[i + 2 * layer_length_] * act_func_first_der_ptr_(z[j][i]);
				dNdwdx[j][i] = act_func_sec_der_ptr_(z[j][i]) * grid_vec[j] * param_vec[i] * param_vec[i + 2 * layer_length_] + act_func_first_der_ptr_(z[j][i]) * param_vec[i + 2 * layer_length_];
				dNdvdx[j][i] = act_func_first_der_ptr_(z[j][i]) * param_vec[i];
				dNdudx[j][i] = act_func_sec_der_ptr_(z[j][i]) * param_vec[i] * param_vec[i + 2 * layer_length_];
			}

		//calculate dNdw, dNdv, dNdu matrices
		for (int j = 0; j < grid_vec.size(); j++) {
			double c = (N[j] + grid_vec[j] * dNdx[j] - f[j]);
			for (int i = 0; i < layer_length_; i++) {			
				dEdw[j][i] = 2 * c * (dNdw[j][i] + grid_vec[j] * dNdwdx[j][i] - dfdu[j] * grid_vec[j] * dNdw[j][i]);
				dEdv[j][i] = 2 * c * (dNdv[j][i] + grid_vec[j] * dNdvdx[j][i] - dfdu[j] * grid_vec[j] * dNdv[j][i]);
				dEdu[j][i] = 2 * c * (dNdu[j][i] + grid_vec[j] * dNdudx[j][i] - dfdu[j] * grid_vec[j] * dNdu[j][i]);
			}
		}

		//---------------------------------------------------------------

		//calculate gradients 
		for (int i = 0; i < layer_length_; i++) {
			for (int j = 0; j < grid_vec.size(); j++) {
				grad_vec[i] += dEdw[j][i];
				grad_vec[i + layer_length_] += dEdv[j][i];
				grad_vec[i + 2* layer_length_] += dEdu[j][i];
			}
		}

		//---------------------------------------------------------------

		//store previous weights
		for (int i = 0; i < 3 * layer_length_; i++) {
			prev_param_vec[i] = param_vec[i];
		}
		double prev_g = g;

		//store previous loss function output
		double prev_L = L;
		
		//---------------------------------------------------------------

		//fragment step until L < prev_L
		bool stop = false;
		while (!stop) {

			//clear data

			for (int j = 0; j < grid_vec.size(); j++)
				for (int i = 0; i < layer_length_; i++) {
					z[j][i] = 0;
				}

			for (int i = 0; i < grid_vec.size(); i++) {

				N[i] = 0;
				dNdx[i] = 0;
				//f[i] = 0;
			}

			//update weights
			for (int i = 0; i < layer_length_; i++) {
				param_vec[i] = param_vec[i] - rate * grad_vec[i];
				param_vec[i + 2 * layer_length_] = param_vec[i + 2 * layer_length_] - rate * grad_vec[i + layer_length_];
				param_vec[i + layer_length_] = param_vec[i + layer_length_] - rate * grad_vec[i + 2 * layer_length_];
			}
		    //g -= rate * dEdg;

			//calculate z, N, dNdX, f

			//calculate z
			for (int j = 0; j < grid_vec.size(); j++)
				for (int i = 0; i < layer_length_; i++) {
					z[j][i] = param_vec[i] * grid_vec[j] + param_vec[i + layer_length_];
				}

			//calculate net outputs
			for (int j = 0; j < grid_vec.size(); j++)
				for (int i = 0; i < layer_length_; i++) {
					N[j] += act_func_ptr_(z[j][i]) * param_vec[i + 2 * layer_length_] + g;
				}

			//calculate dNdx
			for (int j = 0; j < grid_vec.size(); j++)
				for (int i = 0; i < layer_length_; i++) {
					dNdx[j] += param_vec[i + 2 * layer_length_] * param_vec[i] * act_func_first_der_ptr_(z[j][i]);
				}

			//calculate f(x, N)
			for (int i = 0; i < grid_vec.size(); i++) {
				f[i] = func_ptr_(u0 + grid_vec[i] * N[i], grid_vec[i]);
			}

			L = 0;
			//check loss function
			for (int j = 0; j < grid_vec.size(); j++) {
				L += pow(N[j] + grid_vec[j] * dNdx[j] - f[j], 2);
			}
			//std::cout << "----------------" << std::endl;
			//std::cout << "Loss: " << L <<", prev " << prev_L << std::endl;

			if (is_first_iter) {
				stop = true;
			}
			else {
				if (fabs(L) > fabs(prev_L)) {
					rate *= fragment_rate_;
					//std::cout << "Fragmenting rate..." << std::endl;

					//restore weights
					for (int i = 0; i < layer_length_; i++) {
						param_vec[i] = prev_param_vec[i];
						param_vec[i + layer_length_] = prev_param_vec[i + layer_length_];
						param_vec[i + 2*layer_length_] = prev_param_vec[i + 2*layer_length_];
					}
					g = prev_g;
				}
				else {

					grad = 0;
					//check gradient of the loss function
					for (int i = 0; i < layer_length_; i++) {
						grad += param_vec[i] * param_vec[i] + param_vec[i + 2 * layer_length_] * param_vec[i + 2 * layer_length_] + param_vec[i + layer_length_] * param_vec[i + layer_length_];
					}

					grad = sqrt(grad);
					//std::cout << "|grad| =  " << grad << std::endl;

					/*for(int i = 0; i < layer_length_; i++)
						std::cout << "grad(w" << i << ", v" << i << ", u" << i << ") = " << "(" << w_grad[i] << ", " << v_grad[i] << ", " << u_grad[i] << ")" << std::endl;*/
					
					rate *= inv_fragment_rate_;
					stop = true;
				}
			}		
				
		}
		
		bool a = ((fabs(grad) < eps && !is_first_iter) && !(fabs(L - prev_L) < L_threshhold));
		bool b = fabs(L - prev_L) < L_threshhold;
   
		if (a || b) {
			stop_flag = true;
			/*if (a) {
				std::cout << "----------------" << std::endl;
				std::cout << "Gradient precision reached. Finishing..." << std::endl;
			}				
			else
				if (b) {
					std::cout << "----------------" << std::endl;
					std::cout << "Endless loop. Finishing..." << std::endl;
				}	*/				
		}
			
			
		is_first_iter = false;

		it_count++;
	}

	//save new weights
	for (int i = 0; i < layer_length_; i++) {
		net_arr_ptr_[i]->set_weight(1, param_vec[i]);
		net_arr_ptr_[i]->set_bias(param_vec[i + layer_length_]);
		output_neuron_ptr_->set_weight(i+1, param_vec[i + 2 * layer_length_]);
	}

	/*std::cout << "----------------" << std::endl;
	std::cout << "Iterations: " << it_count << std::endl;
	std::cout << "----------------" << std::endl;*/

	loss = L;
}
