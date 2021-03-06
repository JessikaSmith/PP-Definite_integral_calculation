#include <iostream>
#include <stdlib.h>  
#include <omp.h>
#include <math.h> 
#include <chrono>

using namespace std;

static long npoints = 100;


double f(float x) {
	return 1 / pow(x,2) * pow(sin(1 / x),2);
}

void serial_time(double a, double b, double eps, float step) {

	double x_i, dd, J_1, t_1, t_2, dt;
	double edges = (f(a) + f(b)) / 2;
	//double* time = new double[num_of_time_samples];
	double J_2 = 10000000;

	chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now();

	for (int n = 1; n < 10000000000; n > 100000000 ? n += 10000 : n *= step ) {
		dd = (b - a) / n;
		J_1 = J_2;
		J_2 = 0;
		for (int i = 1; i < n-1; i++) {
			x_i = a + dd*i;
			J_2 += f(x_i);
		}
		J_2 += edges;
		J_2 = dd*J_2;
		//cout << log10(abs((J_1 - J_2) / J_2)) <<" "<<J_2<< " and n = " << n << " points." << endl;
		if (log10(abs((J_1 - J_2)/J_2)) <= eps){
			cout << "Result for serial program: " << J_2 << " and n = " << n << " points.";
			break;
		}
	}

	chrono::high_resolution_clock::time_point end = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::seconds>(end - start).count();

	cout << " | Computational time: " << duration << endl;
}


void parallel_time(double a, double b, double eps, float step, int num_of_threads, int num_of_time_samples = 1) {

	double x_i, dd, J_1, t_1, t_2, dt;
	double edges = (f(a) + f(b)) / 2;
	double J_2 = 10000000;
	int i;

	t_1 = omp_get_wtime();

	for (int n = 1; n < 10000000000; n > 100000000 ? n += 100000 : n *= step) {
		dd = (b - a) / n;
		J_1 = J_2;
		J_2 = 0;
		omp_lock_t lock;
		omp_init_lock(&lock);
		#pragma omp parallel for schedule(static) private(x_i) shared(a, dd, n) num_threads(num_of_threads)
		for (i = 1; i < n - 1; i++) {
				x_i = a + dd*i;
				omp_set_lock(&lock);
				J_2 += f(x_i);
				omp_unset_lock(&lock);
		}
		omp_destroy_lock(&lock);
		J_2 += edges;
		J_2 = dd*J_2;
		// cout << log10(abs((J_1 - J_2) / J_2)) <<" "<< J_2 << endl;
		if (log10(abs((J_1 - J_2) / J_2)) <= eps) {
			cout << "Result for " << num_of_threads << " threads: " << J_2 << " and n = " << n << " points.";
			break;
		}
	}
	t_2 = omp_get_wtime();
	dt = t_2 - t_1;

	cout << " | Computational time: " << dt << endl;
}

int main() {

	const int num_of_elem = 7;
	double a[num_of_elem] = { 0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10 };
	double b[num_of_elem] = { 0.0001, 0.001, 0.01, 0.1, 1, 10, 100 };
	double eps[num_of_elem] = { -8, -10, -11, -11, -11, -11, -11 };
	float initial_step = 2;

	// parallel_time(a[4], b[4], eps[4], initial_step, 2);
	for (int i = 6; i > 0; i--) {
		cout << "===========  a = " << a[i] << ", b = " << b[i] << " ===========" << endl;
		serial_time(a[i], b[i], eps[i], initial_step);
		cout << endl;
	}
}