#include "utils.hpp"

void init(int t, double& b) {
	switch (t)
	{
	case 1:
		b = -1;
		break;
	case 2:
		b = -0.24;
		break;
	case 3:
		b = 0;
		break;
	case 4:
		b = -1.4;
		break;
	}
}
