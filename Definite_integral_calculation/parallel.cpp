#include <iostream>
#include <stdlib.h>  
#include <omp.h>
#include <math.h> 

using namespace std;

static long npoints = 100;

double f(float x) {
	return 1.0 / pow(x, 2) * pow(sin(1.0 / x), 2);
}

int main() {
	cout << omp_get_max_threads() << endl;
	double a = 0.00001;
	double b = 0.0001;
	double x_i;
	double eps = 2.77E-7;
	// iteratively getting number of steps for
	// a particular epsilon
	double dd;
	double J_1;
	double J_2 = 10000000;
	double real_J = (double)1/ (double)4*((double)2*((b - a) / (a*b)) + (double)sin((double)2 / b) - (double)sin((double)2/ a));
	cout << real_J;
	double edges = (f(a) + f(b)) / 2;
	//#pragma omp parallel for schedule(static) default(none) \
		shared(a) private()


	for (int n = 3000000; n < 10000000000; n += 1000) {
		dd = (b - a) / n;
		J_1 = J_2;
		J_2 = 0;
//#pragma omp parallel for
		for (int i = 1; i < n - 1; i++) {
			x_i = a + dd*i;
			#pragma omp atomic
			J_2 += f(x_i);
		}
		J_2 += edges;
		J_2 = dd*J_2;
		//if (real_J-J_2)
		if (abs(J_1 - J_2) <= eps*abs(J_2)) break;
	}
	//cout << J_2 << endl
	cout << "Approximated result: " << J_2 << endl;
	cout << "With precision epsilon = " << eps << endl;
}