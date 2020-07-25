#include <iostream>
#include "neural_net_first_order.hpp"
#include "euler.hpp"
#include "runge.hpp"
#include "examples.hpp"
#include "utils.hpp"
#include <chrono>
#include <fstream>
#include <string>
#include <algorithm>
#include <iomanip> 

double eps = 0.0001;

void func(int t, double u_init, double(*f)(double u, double x), double(*f_der)(double u, double x), double(*f_an)(double x), double b ,int param_num, bool analyze, bool analyze_old_meth) {
	if (analyze) {
		std::ofstream analyze;
		std::ofstream plot_layers_rel_err;
		std::ofstream plot_layers_loss;
		analyze.open("analyze.txt");
		plot_layers_rel_err.open("plot_layers_rel_err.txt");
		plot_layers_loss.open("plot_layers_loss.txt");

		analyze << std::fixed << std::setprecision(15);
		plot_layers_rel_err << std::fixed << std::setprecision(15);
		plot_layers_loss << std::fixed << std::setprecision(15);

		double init_step = 0.4;
		for (int k = 1; k <= 6; k++) {

			//grid
			std::vector<double> t_vec;

			double step = init_step * pow(0.5,k);
			int grid_size = (int)(b / step) + 1;

			std::cout << grid_size << " ";

			for (int i = 0; i < grid_size; i++)
				t_vec.push_back(i*step);

			//analyze results w/ different layer length
			NeuralNetFirstOrder net_no_need(2, t);
			
			analyze << std::fixed << std::setprecision(6);
			analyze << "output bias " << net_no_need.get_output_bias() << std::endl;
			analyze << "Grid size: " << grid_size << std::endl;
			analyze << "h = " << step << std::endl;
			analyze << "[" << t_vec[0] << ", " << t_vec[t_vec.size() - 1] << "]" << std::endl;
			analyze << "eps = " << eps << std::endl;
			analyze << "Layer length -- Loss -- Max err -- Iterations -- Parameters num. -- microsec" << std::endl;


			plot_layers_rel_err << std::fixed << std::setprecision(6);
			plot_layers_rel_err << "output bias " << net_no_need.get_output_bias() << std::endl;
			plot_layers_rel_err << "Grid size: " << grid_size << std::endl;
			plot_layers_rel_err << "h = " << step << std::endl;
			plot_layers_rel_err << "[" << t_vec[0] << ", " << t_vec[t_vec.size() - 1] << "]" << std::endl;
			plot_layers_rel_err << "eps = " << eps << std::endl;

			plot_layers_loss << std::fixed << std::setprecision(6);
			plot_layers_loss << "output bias " << net_no_need.get_output_bias() << std::endl;
			plot_layers_loss << "Grid size: " << grid_size << std::endl;
			plot_layers_loss << "h = " << step << std::endl;
			plot_layers_loss << "[" << t_vec[0] << ", " << t_vec[t_vec.size() - 1] << "]" << std::endl;
			plot_layers_loss << "eps = " << eps << std::endl;

			analyze << std::fixed << std::setprecision(15);
			plot_layers_loss << std::fixed << std::setprecision(15);
			plot_layers_rel_err << std::fixed << std::setprecision(15);
			for (int i = 1; i <= 40; i++) {
				int it_count = 0;
				double loss = 0;

				NeuralNetFirstOrder net(i, t);
				auto start = std::chrono::steady_clock::now();
				net.train(u_init, f, f_der, eps, t_vec, it_count, loss);
				auto end = std::chrono::steady_clock::now();

				double time_nano = 0;
				double time_micro = 0;
				double time_milli = 0;
				double time_sec = 0;
				/*for (int d = 0; d < 1000; d++) {
					NeuralNetFirstOrder test_net(i);
					auto start = std::chrono::steady_clock::now();
					test_net.train(u_init, f, f_der, eps, t_vec, it_count, loss);
					auto end = std::chrono::steady_clock::now();

					time_nano += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
					time_micro += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
					time_milli += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
					time_sec += std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
				}

				time_nano /= 1000;
				time_micro /= 1000;
				time_milli /= 1000;
				time_sec /= 1000;*/

				double max_err = 0;
				double max_rel_err = 0;
				for (int j = 0; j < t_vec.size(); j++) {
					double f_net = u_init + t_vec[j] * net.calc_output(t_vec[j]);
					double err = fabs(f_net - f_an(t_vec[j]));
					double rel_err = err / fabs(f_an(t_vec[j]));
					if (j == 0) {
						max_err = err;
						max_rel_err = rel_err;
					}
					else {
						if (max_err < err)
							max_err = err;

						if (max_rel_err < rel_err)
							max_rel_err = rel_err;
					}
						
				}
				
				plot_layers_rel_err << i << " " << max_rel_err << std::endl;
				plot_layers_loss << i << " " << loss << std::endl;
				analyze << i << " | " << loss << " | " << max_err << " | " << it_count << " | " << i * 3  << " | " << time_micro << std::endl;

			}
		
			analyze << std::endl << std::endl << std::endl;
			plot_layers_rel_err << std::endl << std::endl << std::endl;
			plot_layers_loss << std::endl << std::endl << std::endl;

		}
		analyze.close();
		std::ifstream t("analyze.txt");
		std::string str((std::istreambuf_iterator<char>(t)),
			std::istreambuf_iterator<char>());
		std::replace(str.begin(), str.end(), '.', ',');
		t.close();
		analyze.open("analyze.txt");
		analyze << str;	
		analyze.close();


		plot_layers_rel_err.close();
		std::ifstream t2("plot_layers_rel_err.txt");
		std::string str2((std::istreambuf_iterator<char>(t2)),
			std::istreambuf_iterator<char>());
		std::replace(str2.begin(), str2.end(), '.', ',');
		t2.close();
		plot_layers_rel_err.open("plot_layers_rel_err.txt");
		plot_layers_rel_err << str2;
		plot_layers_rel_err.close();


		plot_layers_loss.close();
		std::ifstream t3("plot_layers_rel_err.txt");
		std::string str3((std::istreambuf_iterator<char>(t3)),
			std::istreambuf_iterator<char>());
		std::replace(str3.begin(), str3.end(), '.', ',');
		t3.close();
		plot_layers_loss.open("plot_layers_rel_err.txt");
		plot_layers_loss << str3;
		plot_layers_loss.close();


	}
	

	//grid
	std::vector<double> t_vec;

	double step = 0.1;
	int grid_size = (int)(b / step) + 1;
	std::cout << grid_size << std::endl;

	for (int i = 0; i < grid_size; i++)
		t_vec.push_back(i*step);

	int it_count = 0;
	double loss = 0;

	NeuralNetFirstOrder net(param_num, t);
	net.train(u_init, f, f_der, eps, t_vec, it_count, loss);

	std::vector<double> euler_r;
	std::vector<double> runge_r;
	euler(u_init, f, step, t_vec, euler_r);
	runge_kutta4(u_init, f, step, t_vec, runge_r);


	if (analyze_old_meth) {

		std::ofstream myfile;
		myfile.open("comparison.txt");
		//myfile << std::scientific << std::setprecision(15);
		//myfile << std::fixed << std::setprecision(15);

		double init_step = 0.4;
		for (int k = 1; k <= 6; k++) {

			//grid
			std::vector<double> t_vec;

			double step = init_step * pow(0.5, k);
			int grid_size = (int)(b / step) + 1;

			std::cout << grid_size << " ";

			for (int i = 0; i < grid_size; i++)
				t_vec.push_back(i*step);

			//analyze results w/ different layer length
			NeuralNetFirstOrder net_no_need(2, t);
			myfile << "output bias " << net_no_need.get_output_bias() << std::endl;
			myfile << "Grid size: " << grid_size << std::endl;
			myfile << "h = " << step << std::endl;
			myfile << "[" << t_vec[0] << ", " << t_vec[t_vec.size() - 1] << "]" << std::endl;
			myfile << "eps = " << eps << std::endl;
			myfile << "Layer length -- Loss -- Max err -- Iterations -- Parameters num. -- microsec" << std::endl;

			double time_nano_euler = 0;
			double time_nano_runge = 0;
			double time_nano_net = 0;
			int it_count1 = 0;
			for (int d = 0; d < 1; d++) {
				std::vector<double> euler_res;
				std::vector<double> runge_res;

				NeuralNetFirstOrder net_test(param_num, t);
				auto start0 = std::chrono::steady_clock::now();
				net_test.train(u_init, f, f_der, eps, t_vec, it_count1, loss);
				auto end0 = std::chrono::steady_clock::now();

				auto start1 = std::chrono::steady_clock::now();
				euler(u_init, f, step, t_vec, euler_res);
				auto end1 = std::chrono::steady_clock::now();

				auto start2 = std::chrono::steady_clock::now();
				runge_kutta4(u_init, f, step, t_vec, runge_res);
				auto end2 = std::chrono::steady_clock::now();

				time_nano_euler += std::chrono::duration_cast<std::chrono::nanoseconds>(end1 - start1).count();
				time_nano_runge += std::chrono::duration_cast<std::chrono::nanoseconds>(end2 - start2).count();
				time_nano_net += std::chrono::duration_cast<std::chrono::nanoseconds>(end0 - start0).count();
			}

			time_nano_euler /= 10000;
			time_nano_runge /= 10000;
			time_nano_net /= 10000;

			std::vector<double> euler_res;
			std::vector<double> runge_res;

			euler(u_init, f, step, t_vec, euler_res);
			runge_kutta4(u_init, f, step, t_vec, runge_res);

			double max_err_net = 0;
			double max_rel_err_net = 0;
			for (int j = 0; j < t_vec.size(); j++) {
				double f_net = u_init + t_vec[j] * net.calc_output(t_vec[j]);
				double err = fabs(f_net - f_an(t_vec[j]));
				double rel_err = err / fabs(f_an(t_vec[j]));
				if (j == 0) {
					max_err_net = err;
					max_rel_err_net = rel_err;
				}
				else {
					if (max_err_net < err)
						max_err_net = err;

					if (max_rel_err_net < rel_err)
						max_rel_err_net = rel_err;
				}

			}

			double max_err_euler = 0;
			double max_rel_err_euler = 0;
			for (int j = 0; j < t_vec.size(); j++) {
				double f_net = euler_res[j];
				double err = fabs(f_net - f_an(t_vec[j]));
				double rel_err = err / fabs(f_an(t_vec[j]));
				if (j == 0) {
					max_err_euler = err;
					max_rel_err_euler = rel_err;
				}
				else {
					if (max_err_euler < err)
						max_err_euler = err;

					if (max_rel_err_euler < rel_err)
						max_rel_err_euler = rel_err;
				}

			}

			double max_err_runge = 0.0;
			double max_rel_err_runge = 0.0;
			int max_index = 0;
			for (int j = 0; j < t_vec.size(); j++) {
				double f_net = runge_res[j];
				double err = fabs(f_net - f_an(t_vec[j]));
				double rel_err = err / fabs(f_an(t_vec[j]));
				if (j == 0) {
					max_err_runge = err;
					max_rel_err_runge = rel_err;
					max_index = j;
				}
				else {
					if (max_err_runge < err) {
						max_err_runge = err;
						max_index = j;
					}
						

					if (max_rel_err_net < rel_err) {
						max_rel_err_runge = rel_err;
						//max_index = j;
					}
						
				}
			}

			max_rel_err_runge = max_err_runge / t_vec[max_index];

			myfile << std::fixed;
			myfile << "output bias " << net.get_output_bias() << std::endl;
			myfile << "neurons in layer " << param_num << std::endl;
			myfile << "Grid size: " << grid_size << std::endl;
			myfile << "h = " << step << std::endl;
			myfile << "[" << t_vec[0] << ", " << t_vec[t_vec.size() - 1] << "]" << std::endl << std::endl;

			myfile << std::scientific << std::setprecision(15);
			myfile << "NNet -- Euler -- Runge" << std::endl;
			myfile << "Max abs err.:" << std::endl;
			myfile << max_err_net << " | " << max_err_euler << " | " << max_err_runge << std::endl;
			myfile << "Max rel err.:" << std::endl;
			myfile << max_rel_err_net << " | " << max_rel_err_euler << " | " << max_rel_err_runge << std::endl;
			myfile << "Val: " << t_vec[max_index] << std::endl;

			myfile << std::fixed << std::setprecision(6);
			myfile << "nanosec" << std::endl;
			myfile << time_nano_net << " | " << time_nano_euler << " | " << time_nano_runge << std::endl;
			myfile << "iter: " << it_count1 << std::endl;

			myfile << std::endl << std::endl << std::endl;
		}

		myfile.close();

	}

	std::cout << "----------------" << std::endl;
	std::cout << "Iterations: " << it_count << std::endl;
	std::cout << "Loss: " << loss << std::endl;
	std::cout << "----------------" << std::endl;


	//save results to file
	std::ofstream myfile_coord;
	std::ofstream myfile;
	std::ofstream myfile1;
	std::ofstream myfile2;
	std::ofstream myfile3;
	myfile_coord.open("coord" + std::to_string(t) + ".txt");
	myfile.open("plot" + std::to_string(t) + ".txt");
	myfile1.open("plot_an" + std::to_string(t) + ".txt");
	myfile2.open("plot_euler" + std::to_string(t) + ".txt");
	myfile3.open("plot_runge" + std::to_string(t) + ".txt");

	myfile << std::fixed << std::setprecision(15);
	myfile1 << std::fixed << std::setprecision(15);

	for (int i = 0; i < t_vec.size(); i++) {
		myfile_coord << t_vec[i] << std::endl;
		myfile << u_init + t_vec[i] * net.calc_output(t_vec[i]) << std::endl;
		myfile1 << f_an(t_vec[i]) << std::endl;
		myfile2 << euler_r[i] << std::endl;
		myfile3 << runge_r[i] << std::endl;
	}

	myfile_coord.close();
	myfile.close();
	myfile1.close();
	myfile2.close();
	myfile3.close();
}

int main(){
	std::cout << std::fixed << std::setprecision(15);
	func(1, u0_1, f1, f_der1, f_an1, 1, 4, false, false);
	func(2, u0_2, f2, f_der2, f_an2, 1, 4, false, false);
	func(4, u0_4, f4, f_der4, f_an4, 1, 4, false, false);

	system("pause");
	return 1;
}