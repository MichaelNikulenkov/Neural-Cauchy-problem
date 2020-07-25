#include <iostream>
#include <math.h>
#include <vector>
#include <functional>
#include <fstream>

#define EPS 1e-10
#define CONST_K 20.0
#define CONST_M 0.3
#define CONST_SIGMA 10.0
#define CONST_R 28.0
#define CONST_B 8.0 / 3.0
#define DEFAULT_VALUE 2
#define ADAMS_CONST 4
#define CON_A 0.4
#define CON_B 1.5

using std::vector;
using std::function;
using std::ofstream;

typedef vector<function<double(double, double, double)>> vecOfExpFunc;
typedef vector<function<double(double, double, double, double, double)>> vecOfImpFunc;

vector<double> vecDifference(vector<double> V1, vector<double> V2, double devidedBy) {
	vector<double> result(3, 0);
	for (size_t i = 0; i < V1.size(); ++i) {
		result[i] = (V1[i] - V2[i]) / devidedBy;
	}
	return result;
}

double getNorm(vector<double> V1) {
	double max = EPS;

	for (size_t i = 0; i < V1.size(); ++i) {
		if (fabs(V1[i]) > max) {
			max = fabs(V1[i]);
		}
	}
	return max;
}

void assignFunctions(vecOfExpFunc& explicitF, vecOfImpFunc& implicitF, vecOfImpFunc& schemaF, size_t n) {
	switch (n) {
	case 0: {
		explicitF.emplace_back([](double u1, double u2, double u3) {return u2; });
		explicitF.emplace_back([](double u1, double u2, double u3) {return -u1 * (CONST_K / CONST_M); });
		implicitF.emplace_back([](double u1_Cur, double u1_Prev, double u2_Cur, double blank, double tau) {return u1_Cur - u1_Prev - tau * u2_Cur; });
		implicitF.emplace_back([](double u2_Cur, double u2_Prev, double u1_Cur, double blank, double tau) {return u2_Cur - u2_Prev + tau * u1_Cur * (CONST_K / CONST_M); });
		schemaF.emplace_back([](double u1_Cur, double u1_Prev, double u2_Cur, double u2_Prev, double tau) {return u1_Cur - u1_Prev - (tau / 2) * (u2_Cur + u2_Prev); });
		schemaF.emplace_back([](double u2_Cur, double u2_Prev, double u1_Cur, double u1_Prev, double tau) {return u2_Cur - u2_Prev + (tau / 2) * (CONST_K / CONST_M) * (u1_Cur + u1_Prev); });
		break;
	}
	case 1: {
		explicitF.emplace_back([](double x, double y, double z) {return 2.0 * x + y * y - 1.0; });
		explicitF.emplace_back([](double x, double y, double z) {return 6.0 * x - y * y + 1.0; });
		implicitF.emplace_back([](double x_Cur, double x_Prev, double y_Cur, double blank, double tau) {return x_Cur - x_Prev - tau * (2.0 * x_Cur + y_Cur * y_Cur - 1.0); });
		implicitF.emplace_back([](double y_Cur, double y_Prev, double x_Cur, double blank, double tau) {return y_Cur - y_Prev - tau * (6.0 * x_Cur - y_Cur * y_Cur + 1.0); });
		schemaF.emplace_back([](double x_Cur, double x_Prev, double y_Cur, double y_Prev, double tau) {return x_Cur - x_Prev - (tau / 2) * (2.0 * x_Cur + y_Cur * y_Cur - 1.0 + 2.0 * x_Prev + y_Prev * y_Prev - 1.0); });
		schemaF.emplace_back([](double y_Cur, double y_Prev, double x_Cur, double x_Prev, double tau) {return y_Cur - y_Prev - (tau / 2) * (6.0 * x_Cur - y_Cur * y_Cur + 1.0 + 6.0 * x_Prev - y_Prev * y_Prev + 1.0); });
		break;
	}
	case 2: {
		explicitF.emplace_back([](double x, double y, double z) { return 1.0 - x * x - y * y; });
		explicitF.emplace_back([](double x, double y, double z) { return 2.0 * x; });
		implicitF.emplace_back([](double x_Cur, double x_Prev, double y_Cur, double blank, double tau) {return x_Cur - x_Prev - tau * (1.0 - x_Cur * x_Cur - y_Cur * y_Cur); });
		implicitF.emplace_back([](double y_Cur, double y_Prev, double x_Cur, double blank, double tau) {return y_Cur - y_Prev - tau * (2.0 * x_Cur); });
		schemaF.emplace_back([](double x_Cur, double x_Prev, double y_Cur, double y_Prev, double tau) {return x_Cur - x_Prev - (tau / 2) * (1.0 - x_Cur * x_Cur - y_Cur * y_Cur + 1.0 - x_Prev * x_Prev - y_Prev * y_Prev); });
		schemaF.emplace_back([](double y_Cur, double y_Prev, double x_Cur, double x_Prev, double tau) {return y_Cur - y_Prev - (tau / 2) * (2.0 * x_Cur + 2.0 * x_Prev); });
		break;
	}
	case 3: {
		explicitF.emplace_back([](double x, double y, double z) {return CONST_SIGMA * (y - x); });
		explicitF.emplace_back([](double x, double y, double z) {return x * (CONST_R - z) - y; });
		explicitF.emplace_back([](double x, double y, double z) {return x * y - CONST_B * z; });
		//implicitF.emplace_back( [](double x_Cur, double x_Prev, double y_Cur, double z_Cur, double tau) {return x_Cur - x_Prev - tau * CONST_SIGMA *(y_Cur - x_Cur);} );
		//implicitF.emplace_back( [](double y_Cur, double y_Prev, double x_Cur, double z_Cur, double tau) {return y_Cur - y_Prev - tau * (x_Cur * (CONST_R - z_Cur) - y_Cur);} );
		//implicitF.emplace_back( [](double z_Cur, double z_Prev, double x_Cur, double y_Cur, double tau) {return z_Cur - z_Prev - tau * (x_Cur * y_Cur - CONST_B * z_Cur);} );
		break;
	}
	case 4: {
		explicitF.emplace_back([](double x, double y, double z) {return CON_A - (CON_B + 1) * x + (x * x * y); });
		explicitF.emplace_back([](double x, double y, double z) {return CON_B * x - (x * x * y); });
		implicitF.emplace_back([](double x_Cur, double x_Prev, double y_Cur, double blank, double tau) {return x_Cur - x_Prev - tau * (CON_A - (CON_B + 1) * x_Cur + (x_Cur * x_Cur * y_Cur)); });
		implicitF.emplace_back([](double y_Cur, double y_Prev, double x_Cur, double blank, double tau) {return y_Cur - y_Prev - tau * (CON_B * x_Cur - (x_Cur * x_Cur * y_Cur)); });
		schemaF.emplace_back([](double x_Cur, double x_Prev, double y_Cur, double y_Prev, double tau) {return x_Cur - x_Prev - (tau / 2) * (CON_A - (CON_B + 1) * x_Cur + (x_Cur * x_Cur * y_Cur) + CON_A - (CON_B + 1) * x_Prev + (x_Prev * x_Prev * y_Prev)); });
		schemaF.emplace_back([](double y_Cur, double y_Prev, double x_Cur, double x_Prev, double tau) {return y_Cur - y_Prev - (tau / 2) * (CON_B * x_Cur - (x_Cur * x_Cur * y_Cur) + CON_B * x_Prev - (x_Prev * x_Prev * y_Prev)); });
	}
	default: {return; }
	}
}

vector<double> Gauss(vector<vector<double>> A, vector<double> b, size_t n) {
	vector<vector<double>> A_Tmp(n, vector<double>(n + 1, 0));

	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < n; ++j) {
			A_Tmp[i][j] = A[i][j];
		}
	}

	for (size_t i = 0; i < n; ++i) {
		A_Tmp[i][n] = b[i];
	}

	double max = 0, value_Tmp = 0;
	size_t jMax = 0;
	bool flag1 = false, flag2 = true;

	for (size_t i = 0; i < n; ++i) {
		if (flag2) {
			max = fabs(A_Tmp[i][i]);

			for (size_t j = i + 1; j < n; ++j) {
				flag1 = false;

				if (fabs(A_Tmp[j][i]) > fabs(max)) {
					max = fabs(A_Tmp[j][i]);
					jMax = j;
					flag1 = true;
				}

				if (flag1) {
					for (size_t k = 0; k < n + 1; ++k) {
						std::swap(A_Tmp[i][k], A_Tmp[jMax][k]);
					}
				}
			}

			for (size_t j = i + 1; j < n; ++j) {
				value_Tmp = A_Tmp[j][i] / A_Tmp[i][i];

				for (size_t k = 0; k < n; ++k) {
					A_Tmp[j][k] = A_Tmp[j][k] - value_Tmp * A_Tmp[i][k];
				}
			}
		}

		if (fabs(A_Tmp[i][i]) < EPS) {
			flag2 = false;
		}
	}

	vector<double> x(n, 0);

	if (flag2) {
		for (int i = n - 1; i >= 0; --i) {

			double sum = 0;

			for (size_t j = n - 1; j > i; --j) {
				sum += x[j] * A_Tmp[i][j];
			}

			x[i] = (A_Tmp[i][n] - sum) / A_Tmp[i][i];

		}
	}
	else {
		for (size_t i = 0; i < n; ++i) {
			x[i] = 0;
		}
	}
	return x;
}

vector<vector<double>> Jacobi(double x, double x_Prev, double y, double y_Prev, double z, double z_Prev, double tau, size_t n, vecOfImpFunc implicitF) {
	switch (n) {
	case 0:
	case 1:
	case 2: {
		vector<vector<double>> F(3, vector<double>(2, 0));
		F[0][0] = (implicitF[0](x + EPS, x_Prev, y, y_Prev, tau) - implicitF[0](x, x_Prev, y, y_Prev, tau)) / EPS;
		F[0][1] = (implicitF[0](x, x_Prev, y + EPS, y_Prev, tau) - implicitF[0](x, x_Prev, y, y_Prev, tau)) / EPS;
		F[1][0] = (implicitF[1](y + EPS, y_Prev, x, x_Prev, tau) - implicitF[1](y, y_Prev, x, x_Prev, tau)) / EPS;
		F[1][1] = (implicitF[1](y, y_Prev, x + EPS, x_Prev, tau) - implicitF[1](y, y_Prev, x, x_Prev, tau)) / EPS;
		return F;
	}
	case 3: {
		vector<vector<double>> F(3, vector<double>(3, 0));
		F[0][0] = (implicitF[0](x + EPS, x_Prev, y, z, tau) - implicitF[0](x, x_Prev, y, z, tau)) / EPS;
		F[0][1] = (implicitF[0](x, x_Prev, y + EPS, z, tau) - implicitF[0](x, x_Prev, y, z, tau)) / EPS;
		F[0][2] = (implicitF[0](x, x_Prev, y, z + EPS, tau) - implicitF[0](x, x_Prev, y, z, tau)) / EPS;
		F[1][0] = (implicitF[1](x + EPS, x_Prev, y, z, tau) - implicitF[1](x, x_Prev, y, z, tau)) / EPS;
		F[1][1] = (implicitF[1](x, x_Prev, y + EPS, z, tau) - implicitF[1](x, x_Prev, y, z, tau)) / EPS;
		F[1][2] = (implicitF[1](x, x_Prev, y, z + EPS, tau) - implicitF[1](x, x_Prev, y, z, tau)) / EPS;
		F[2][0] = (implicitF[2](x + EPS, x_Prev, y, z, tau) - implicitF[2](x, x_Prev, y, z, tau)) / EPS;
		F[2][1] = (implicitF[2](x, x_Prev, y + EPS, z, tau) - implicitF[2](x, x_Prev, y, z, tau)) / EPS;
		F[2][2] = (implicitF[2](x, x_Prev, y, z + EPS, tau) - implicitF[2](x, x_Prev, y, z, tau)) / EPS;
		return F;
	}
	default: return vector<vector<double>>(2, vector<double>(2, 0));
	}
}

vector<double> sysNewton(double x_Prev, double y_Prev, double z_Prev, vecOfImpFunc implicitF, double tau, size_t n) {
	vector<vector<double>> jacobiMatrix(3, vector<double>(3, 0));
	double x_Tmp = 0, y_Tmp = 0, z_Tmp = 0, k = 0, x = x_Prev, y = y_Prev, z = z_Prev;
	size_t iter = 0, iter_Tmp = 0, iter_Diff = 0;
	bool notStop = true;
	vector<double> nextVector(3, 0);

	while (notStop) {
		vector<double> Y(3, 0);
		jacobiMatrix = Jacobi(x, x_Prev, y, y_Prev, z, z_Prev, tau, n, implicitF);

		switch (n) {
		case 0:
		case 1:
		case 4:
		case 2: { Y[0] = -implicitF[0](x, x_Prev, y, y_Prev, tau); Y[1] = -implicitF[1](y, y_Prev, x, x_Prev, tau); Y[2] = 0; break; }
		case 3: { Y[0] = -implicitF[0](x, x_Prev, y, z, tau); Y[1] = -implicitF[1](y, y_Prev, x, z, tau); Y[2] = -implicitF[2](z, z_Prev, x, y, tau); break; }
		default: {return vector<double>(3, 0); }
		}

		n == 3 ? Y = Gauss(jacobiMatrix, Y, 3) : Gauss(jacobiMatrix, Y, 2);

		++iter;
		x_Tmp = x; y_Tmp = y; z_Tmp = z; x += Y[0]; y += Y[1]; z += Y[2];
		iter_Diff = iter - iter_Tmp;

		if (sqrt((x - x_Tmp) * (x - x_Tmp) + (y - y_Tmp) * (y - y_Tmp)) + (z - z_Tmp) * (z - z_Tmp) < EPS || iter_Diff > 30) {
			notStop = false;
			nextVector[0] = x; nextVector[1] = y; nextVector[2] = z;
		}
	}
	return nextVector;
}

void explicitEulerMethod(vecOfExpFunc explicitF, vector<double> initVector, size_t n, double a, double b, double tau) {
	vector<vector<double>> funcY(3, vector<double>(2, 0));
	ofstream file("explicitEulerMethod.dat");

	for (size_t i = 0; i < explicitF.size(); ++i) {
		funcY[i][0] = initVector[i];
		file << funcY[i][0] << "\t\t";
	}

	file << std::endl;
	size_t k = 1;

	while (a + k * tau <= b) {
		for (size_t i = 0; i < explicitF.size(); ++i) {
			funcY[i][1] = funcY[i][0] + tau * explicitF[i](funcY[0][0], funcY[1][0], funcY[2][0]);
			file << funcY[i][1] << "\t\t";
		}

		file << std::endl;

		for (size_t i = 0; i < explicitF.size(); ++i) {
			funcY[i][0] = funcY[i][1];
		}

		++k;
	}
	file.close();
}

void implicitEulerMethod(vecOfImpFunc implicitF, vector<double> initVector, size_t n, double a, double b, double tau) {
	vector<vector<double>> funcY(3, vector<double>(2, 0));
	ofstream file("implicitEulerMethod.dat");

	for (size_t i = 0; i < implicitF.size(); ++i) {
		funcY[i][0] = initVector[i];
		file << funcY[i][0] << "\t\t";
	}

	file << std::endl;
	size_t k = 1;

	while (a + k * tau <= b) {
		vector<double> tmpVector = sysNewton(funcY[0][0], funcY[1][0], funcY[2][0], implicitF, tau, n);

		for (size_t i = 0; i < implicitF.size(); ++i) {
			funcY[i][1] = tmpVector[i];
			file << funcY[i][1] << "\t\t";
		}

		file << std::endl;

		for (size_t i = 0; i < implicitF.size(); ++i) {
			funcY[i][0] = funcY[i][1];
		}

		++k;
	}
	file.close();
}

void symSchemaMethod(vecOfImpFunc schemaF, vector<double> initVector, size_t n, double a, double b, double tau) {
	vector<vector<double>> funcY(3, vector<double>(2, 0));
	ofstream file("symSchemaMethod.dat");

	for (size_t i = 0; i < schemaF.size(); ++i) {
		funcY[i][0] = initVector[i];
		file << funcY[i][0] << "\t\t";
	}

	file << std::endl;
	size_t k = 1;

	while (a + k * tau <= b) {
		vector<double> tmpVector = sysNewton(funcY[0][0], funcY[1][0], funcY[2][0], schemaF, tau, n);

		for (size_t i = 0; i < schemaF.size(); ++i) {
			funcY[i][1] = tmpVector[i];
			file << funcY[i][1] << "\t\t";
		}

		file << std::endl;

		for (size_t i = 0; i < schemaF.size(); ++i) {
			funcY[i][0] = funcY[i][1];
		}

		++k;
	}
	file.close();
}

void rkSecondOrderMethod(vecOfExpFunc explicitF, vector<double> initVector, size_t n, double a, double b, double tau) {
	vector<vector<double>> funcY(3, vector<double>(2, 0));
	ofstream file("rkSecondOrderMethod.dat");

	for (size_t i = 0; i < explicitF.size(); ++i) {
		funcY[i][0] = initVector[i];
		file << funcY[i][0] << "\t\t";
	}

	file << std::endl;
	size_t k = 1;

	while (a + k * tau <= b) {
		vector<double> k1(3, 0);
		vector<double> k2(3, 0);

		for (size_t j = 0; j < explicitF.size(); ++j) {
			k1[j] = (explicitF[j](funcY[0][0], funcY[1][0], funcY[2][0]));
		}

		for (size_t j = 0; j < explicitF.size(); ++j) {
			size_t ind = 0;
			k2[j] = (explicitF[j](funcY[0][0] + tau * k1[ind], funcY[1][0] + tau * k1[ind + 1], funcY[2][0] + tau * k1[ind + 2]));
		}

		for (size_t i = 0; i < explicitF.size(); ++i) {
			funcY[i][1] = funcY[i][0] + tau * (0.5 * k1[i] + 0.5 * k2[i]);
			file << funcY[i][1] << "\t\t";
		}

		file << std::endl;

		for (size_t i = 0; i < explicitF.size(); ++i) {
			funcY[i][0] = funcY[i][1];
		}

		++k;
	}
	file.close();
}

vector<vector<double>> rkFourthOrderMethod(vecOfExpFunc explicitF, vector<double> initVector, size_t n, double a, double b, double tau, bool isAdams = false) {

	ofstream file("rkFourthOrderMethod.dat");
	vector<vector<double>> funcY(3, vector<double>(2, 0));
	for (size_t i = 0; i < explicitF.size(); ++i) {
		funcY[i][0] = initVector[i];
		file << funcY[i][0] << "\t\t";
	}

	vector<vector<double>> adamsVector(3, vector<double>(5, 0));

	for (size_t i = 0; i < explicitF.size(); ++i) {
		adamsVector[i][0] = funcY[i][0];
	}

	file << std::endl;
	size_t k = 1, index = 1;

	while (a + k * tau <= b) {
		vector<double> k1(3, 0);
		vector<double> k2(3, 0);
		vector<double> k3(3, 0);
		vector<double> k4(3, 0);

		for (size_t j = 0; j < explicitF.size(); ++j) {
			k1[j] = (explicitF[j](funcY[0][0], funcY[1][0], funcY[2][0]));
		}

		for (size_t j = 0; j < explicitF.size(); ++j) {
			size_t ind = 0;
			k2[j] = (explicitF[j](funcY[0][0] + (tau / 2) * k1[ind], funcY[1][0] + (tau / 2) * k1[ind + 1],
				funcY[2][0] + (tau / 2) * k1[ind + 2]));
		}

		for (size_t j = 0; j < explicitF.size(); ++j) {
			size_t ind = 0;
			k3[j] = (explicitF[j](funcY[0][0] + (tau / 2) * k2[ind], funcY[1][0] + (tau / 2) * k2[ind + 1],
				funcY[2][0] + (tau / 2) * k2[ind + 2]));
		}

		for (size_t j = 0; j < explicitF.size(); ++j) {
			size_t ind = 0;
			k4[j] = (explicitF[j](funcY[0][0] + tau * k3[ind], funcY[1][0] + tau * k3[ind + 1],
				funcY[2][0] + tau * k3[ind + 2]));
		}

		for (size_t i = 0; i < explicitF.size(); ++i) {
			funcY[i][1] = funcY[i][0] + (tau / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
			file << funcY[i][1] << "\t\t";
			if (isAdams) {
				if (index < 4) {
					adamsVector[i][index] = funcY[i][1];
				}
				else {
					return adamsVector;
				}
			}
		}

		if (isAdams) {
			++index;
		}

		file << std::endl;

		for (size_t i = 0; i < explicitF.size(); ++i) {
			funcY[i][0] = funcY[i][1];
		}

		++k;
	}

	file.close();
	return vector<vector<double>>(3, vector<double>(2, 0));
}

void adamsMethod(vecOfExpFunc explicitF, vector<double> initVector, size_t n, double a, double b, double tau) {
	ofstream file("adamsMethod.dat");

	vector<vector<double>> funcY = rkFourthOrderMethod(explicitF, initVector, n, a, b, tau, true);

	for (size_t i = 0; i < ADAMS_CONST; ++i) {
		for (size_t j = 0; j < explicitF.size(); ++j) {
			file << funcY[j][i] << "\t\t";
		}
		file << std::endl;
	}

	size_t k = 1;

	while (a + k * tau <= b) {
		for (size_t i = 0; i < explicitF.size(); ++i) {
			funcY[i][ADAMS_CONST] = funcY[i][ADAMS_CONST - 1] + (tau / 24) *
				(
					55 * explicitF[i](funcY[0][ADAMS_CONST - 1], funcY[1][ADAMS_CONST - 1], funcY[2][ADAMS_CONST - 1]) -
					59 * explicitF[i](funcY[0][ADAMS_CONST - 2], funcY[1][ADAMS_CONST - 2], funcY[2][ADAMS_CONST - 2]) +
					37 * explicitF[i](funcY[0][ADAMS_CONST - 3], funcY[1][ADAMS_CONST - 3], funcY[2][ADAMS_CONST - 3]) -
					9 * explicitF[i](funcY[0][ADAMS_CONST - 4], funcY[1][ADAMS_CONST - 4], funcY[2][ADAMS_CONST - 4])
					);

			file << funcY[i][ADAMS_CONST] << "\t\t";
		}

		file << std::endl;

		for (size_t i = 0; i < explicitF.size(); ++i) {
			for (size_t index = 0; index < ADAMS_CONST; ++index) {
				funcY[i][index] = funcY[i][index + 1];
			}
		}

		++k;
	}
	file.close();
}

void pkMethod(vecOfExpFunc explicitF, vector<double> initVector, size_t n, double a, double b, double tau) {
	ofstream file("pkMethod.dat");

	vector<vector<double>> funcY = rkFourthOrderMethod(explicitF, initVector, n, a, b, tau, true);

	for (size_t i = 0; i < ADAMS_CONST; ++i) {
		for (size_t j = 0; j < explicitF.size(); ++j) {
			file << funcY[j][i] << "\t\t";
		}
		file << std::endl;
	}

	size_t k = 1;

	while (a + k * tau <= b) {
		vector<double> funcZero(3, 0);
		for (size_t i = 0; i < explicitF.size(); ++i) {
			funcY[i][ADAMS_CONST] = funcY[i][ADAMS_CONST - 1] + (tau / 24) *
				(
					55 * explicitF[i](funcY[0][ADAMS_CONST - 1], funcY[1][ADAMS_CONST - 1], funcY[2][ADAMS_CONST - 1]) -
					59 * explicitF[i](funcY[0][ADAMS_CONST - 2], funcY[1][ADAMS_CONST - 2], funcY[2][ADAMS_CONST - 2]) +
					37 * explicitF[i](funcY[0][ADAMS_CONST - 3], funcY[1][ADAMS_CONST - 3], funcY[2][ADAMS_CONST - 3]) -
					9 * explicitF[i](funcY[0][ADAMS_CONST - 4], funcY[1][ADAMS_CONST - 4], funcY[2][ADAMS_CONST - 4])
					);

		}

		for (size_t i = 0; i < explicitF.size(); ++i) {
			funcZero[i] = explicitF[i](funcY[0][ADAMS_CONST], funcY[1][ADAMS_CONST], funcY[2][ADAMS_CONST]);
		}

		for (size_t i = 0; i < explicitF.size(); ++i) {
			funcY[i][ADAMS_CONST] = funcY[i][ADAMS_CONST - 1] + (tau / 24) *
				(
					9 * funcZero[i] +
					19 * explicitF[i](funcY[0][ADAMS_CONST - 1], funcY[1][ADAMS_CONST - 1], funcY[2][ADAMS_CONST - 1]) -
					5 * explicitF[i](funcY[0][ADAMS_CONST - 2], funcY[1][ADAMS_CONST - 2], funcY[2][ADAMS_CONST - 2]) +
					1 * explicitF[i](funcY[0][ADAMS_CONST - 3], funcY[1][ADAMS_CONST - 3], funcY[2][ADAMS_CONST - 3])
					);

			file << funcY[i][ADAMS_CONST] << "\t\t";
		}

		file << std::endl;

		for (size_t i = 0; i < explicitF.size(); ++i) {
			for (size_t index = 0; index < ADAMS_CONST; ++index) {
				funcY[i][index] = funcY[i][index + 1];
			}
		}

		++k;
	}
	file.close();
}

void rungeRule(vecOfExpFunc explicitF, vector<double> initVector, size_t n, double a, double b, double tau, bool isAdams = false) {

	ofstream file("rungeRule.dat");
	vector<vector<double>> funcY(3, vector<double>(2, 0));

	for (size_t i = 0; i < explicitF.size(); ++i) {
		funcY[i][0] = initVector[i];
		file << funcY[i][0] << "\t\t";
	}

	file << std::endl;
	size_t k = 1;

	while (a + k * tau <= b) {
		vector<double> fullTauValueVector(3, 0);
		vector<double> halfTauValueVector(3, 0);
		do {
			vector<double> k1(3, 0), k2(3, 0), k3(3, 0), k4(3, 0);
			vector<double> _k2(3, 0), _k3(3, 0), _k4(3, 0);

			for (size_t j = 0; j < explicitF.size(); ++j) {
				k1[j] = (explicitF[j](funcY[0][0], funcY[1][0], funcY[2][0]));
			}

			for (size_t j = 0; j < explicitF.size(); ++j) {
				size_t ind = 0;
				k2[j] = (explicitF[j](funcY[0][0] + (tau / 2) * k1[ind], funcY[1][0] + (tau / 2) * k1[ind + 1],
					funcY[2][0] + (tau / 2) * k1[ind + 2]));
				_k2[j] = (explicitF[j](funcY[0][0] + (tau / 4) * k1[ind], funcY[1][0] + (tau / 4) * k1[ind + 1],
					funcY[2][0] + (tau / 4) * k1[ind + 2]));
			}

			for (size_t j = 0; j < explicitF.size(); ++j) {
				size_t ind = 0;
				k3[j] = (explicitF[j](funcY[0][0] + (tau / 2) * k2[ind], funcY[1][0] + (tau / 2) * k2[ind + 1],
					funcY[2][0] + (tau / 2) * k2[ind + 2]));
				_k3[j] = (explicitF[j](funcY[0][0] + (tau / 4) * _k2[ind], funcY[1][0] + (tau / 4) * _k2[ind + 1],
					funcY[2][0] + (tau / 4) * _k2[ind + 2]));
			}

			for (size_t j = 0; j < explicitF.size(); ++j) {
				size_t ind = 0;
				k4[j] = (explicitF[j](funcY[0][0] + tau * k3[ind], funcY[1][0] + tau * k3[ind + 1],
					funcY[2][0] + tau * k3[ind + 2]));
				_k4[j] = (explicitF[j](funcY[0][0] + (tau / 2) * _k3[ind], funcY[1][0] + (tau / 2) * _k3[ind + 1],
					funcY[2][0] + (tau / 2) * _k3[ind + 2]));
			}



			for (size_t i = 0; i < explicitF.size(); ++i) {
				fullTauValueVector[i] = funcY[i][0] + (tau / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
				halfTauValueVector[i] = funcY[i][0] + (tau / 12) * (k1[i] + 2 * _k2[i] + 2 * _k3[i] + _k4[i]);
				//funcY[i][1] = funcY[i][0] + (tau / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
				//file << funcY[i][1] << "\t\t";
			}

			double test = getNorm(vecDifference(halfTauValueVector, fullTauValueVector, 15.0));

			if (test <= 0.0001) {
				for (size_t i = 0; i < explicitF.size(); ++i) {
					funcY[i][1] = halfTauValueVector[i];
					file << funcY[i][1] << "\t\t";
				}
				tau *= 2;
				break;
			}
			else {
				tau /= 2;
			}

		} while (getNorm(vecDifference(halfTauValueVector, fullTauValueVector, 15.0)) > EPS); // getNorm function,


		file << std::endl;

		for (size_t i = 0; i < explicitF.size(); ++i) {
			funcY[i][0] = funcY[i][1];
		}

		++k;
	}

	file.close();
}

//int main() {
//	std::cout << "Введите номер теста\n\n";
//	size_t n = 0;
//	std::cin >> n;
//
//	vecOfExpFunc explicitFunctions;
//	vecOfImpFunc implicitFunctions, schemaFunctions;
//
//	assignFunctions(explicitFunctions, implicitFunctions, schemaFunctions, n);
//
//	std::cout << "Введите границы отрезка\n\n";
//	double a = 0, b = 0, tau = 0.01;
//	std::cin >> a >> b;
//
//	std::cout << "Введите начальные условия\n\n";
//	vector<double> initialVector;
//
//	for (size_t i = 0; i < DEFAULT_VALUE; ++i) {
//		double value = 0;
//		std::cin >> value;
//		initialVector.emplace_back(value);
//	}
//
//	if (n == 3) {
//		double value = 0;
//		std::cin >> value;
//		initialVector.emplace_back(value);
//	}
//
//	//    std::cout << "\n\nВыполняется явный метод Эйлера...\n";
//	//    explicitEulerMethod(explicitFunctions, initialVector, n, a, b, tau);
//	//    std::cout << "Готово!\n\n";
//	//
//	//    std::cout << "Выполняется неявный метод Эйлера...\n";
//	//    implicitEulerMethod(implicitFunctions, initialVector, n, a, b, tau);
//	//    std::cout << "Готово!\n\n";
//	//
//	//    std::cout << "Выполняется симметричная схема...\n";
//	//    symSchemaMethod(schemaFunctions, initialVector, n, a, b, tau);
//	//    std::cout << "Готово!\n\n";
//	//
//	//    std::cout << "Выполняется метод Рунге-Кутты 2 порядка...\n";
//	//    rkSecondOrderMethod(explicitFunctions, initialVector, n, a, b, tau);
//	//    std::cout << "Готово!\n\n";
//	//
//	//    std::cout << "Выполняется метод Адамса-Башфорта 4 порядка...\n";
//	//    adamsMethod(explicitFunctions, initialVector, n, a, b, tau);
//	//    std::cout << "Готово!\n\n";
//	//
//	//    std::cout << "Выполняется метод Прогноза-Коррекции...\n";
//	//    pkMethod(explicitFunctions, initialVector, n, a, b, tau);
//	//    std::cout << "Готово!\n\n";
//
//	//    std::cout << "Выполняется метод Рунге-Кутты 4 порядка...\n";
//	//    rkFourthOrderMethod(explicitFunctions, initialVector, n, a, b, tau);
//	//    std::cout << "Готово!\n\n";
//
//	/*std::cout << "Выполняется правило Рунге...\n";
//	rungeRule(explicitFunctions, initialVector, n, a, b, tau);
//	std::cout << "Готово!\n\n";*/
//
//	return 0;
//}