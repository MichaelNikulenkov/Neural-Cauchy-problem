#pragma once
#include <vector>

class SquareMatrix {
public:
	SquareMatrix(int size, bool is_identity = false);	
	~SquareMatrix();
	void mult_right_vec(const std::vector<double> &vec, std::vector<double> &result);

	int get_size();

	//i - row, j - column
	void set(int i, int j, double val);
	double get(int i, int j);

private:
	int size_ = 0;
	double** matrix_ = nullptr;
};

