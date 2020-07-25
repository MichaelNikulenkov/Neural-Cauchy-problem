#include "square_matrix.hpp"

SquareMatrix::SquareMatrix(int size, bool is_identity) {

	matrix_ = new double*[size];
	for (int i = 0; i < size; i++)
		matrix_[i] = new double[size];
	  
	if (is_identity) {
		for (int i = 0; i < size; i++)
			for (int j = 0; j < size; j++)
				if(i == j)
					matrix_[i][j] = 1;
				else
				    matrix_[i][j] = 0;
	}
	else {
		for (int i = 0; i < size; i++)
			for (int j = 0; j < size; j++)
				matrix_[i][j] = 0;
	}

}

SquareMatrix::~SquareMatrix() {
	for (int i = 0; i < size_; i++)
		delete[] matrix_[i];
	delete[] matrix_;
}

int SquareMatrix::get_size() {
	return size_;
}

double SquareMatrix::get(int i, int j) {
	return matrix_[i][j];
}

void SquareMatrix::set(int i, int j, double val) {
	matrix_[i][j] = val;
}

void SquareMatrix::mult_right_vec(const std::vector<double> &vec, std::vector<double> &result) {
	for (int i = 0; i < size_; i++)
		for (int j = 0; j < size_; j++)
			result[i] += matrix_[i][j] * vec[j];
}