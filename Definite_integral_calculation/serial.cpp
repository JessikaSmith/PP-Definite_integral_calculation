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

void serial_time(double a, double b, double eps, float step, int num_of_time_samples = 1) {

	double x_i, dd, J_1, t_1, t_2, dt;

	//double real_J = (double)1 / (double)4 * ((double)2 * ((b - a) / (a*b)) + (double)sin((double)2 / b) - (double)sin((double)2 / a));
	double edges = (f(a) + f(b)) / 2;
	double* time = new double[num_of_time_samples];
	double J_2;
	for (int t = 1; t < num_of_time_samples + 1; t++) {
		J_2 = 10000000;

		std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

		for (int n = 3000; n < 10000000000; n += step) {
			dd = (b - a) / n;
			J_1 = J_2;
			J_2 = 0;
			for (int i = 1; i < n-1; i++) {
				x_i = a + dd*i;
				J_2 += f(x_i);
			}
			J_2 += edges;
			J_2 = dd*J_2;
			//cout << log10(abs(J_1 - J_2)) << endl;
			//cout << n << endl;
			if (log10(abs(J_1 - J_2)) <= eps) break;
		}

		std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
		time[t] = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
	}
	double minim = time[0];
	for (int y = 1; y < num_of_time_samples; y++) {
		if (time[y] < minim) minim = time[y];
	}
	cout << "Result for serial program: " << J_2 << " | Computational time: " << minim/1000 << endl;

}

void parallel_time(double a, double b, double eps, float step, int num_of_threads, int num_of_time_samples = 1) {

	double x_i, dd, J_1, t_1, t_2, dt;
	//double real_J = (double)1 / (double)4 * ((double)2 * ((b - a) / (a*b)) + (double)sin((double)2 / b) - (double)sin((double)2 / a));
	double edges = (f(a) + f(b)) / 2;
	double* time = new double[num_of_time_samples];
	double J_2;
	int n, i;
	for (int t = 1; t < num_of_time_samples + 1; t++) {
		J_2 = 10000000;

		omp_lock_t lock;
		omp_init_lock(&lock);

		t_1 = omp_get_wtime();
		dd = 0;
		#pragma omp parallel for schedule(static) private(i,n) \
			num_threads(num_of_threads) reduction(+:edges)
		for (n = 3000; n < 10000000000; n += step) {
			#pragma omp atomic
			dd += (b - a) / n;
			J_1 = J_2;
			J_2 = 0;
			for (i = 1; i < n - 1; i++) {
				#pragma omp critical
				x_i = a + dd*i;
				omp_set_lock(&lock);
				J_2 += f(x_i);
				omp_unset_lock(&lock);

			}
			J_2 += edges;
			#pragma omp critical
			J_2 = dd*J_2;
			if (log10(abs(J_1 - J_2)) <= eps) break;
		}

		t_2 = omp_get_wtime();
		dt = t_2 - t_1;
		time[t] = dt;
	}

	double minim = time[0];
	for (int y = 1; y < num_of_time_samples; y++) {
		if (time[y] < minim) minim = time[y];
	}
	cout << "Result for " << num_of_threads << " threads: " << J_2 << " | Computational time: " << minim << endl;

}


int main() {

	const int num_of_elem = 7;
	double a[num_of_elem] = { 0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10 };
	double b[num_of_elem] = { 0.0001, 0.001, 0.01, 0.1, 1, 10, 100 };
	//double eps[num_of_elem] = { -2.77E-11, 1.09E-10, 2.05E-11, -2.22E-12, 8.67E-11, -6.00E-11, -6.00E-11 };
	double eps[num_of_elem] = { -11, -10, -11, -12, -11, -11, -11 };
	float initial_step = 10000;
	int num_of_time_samples = 1;

	for (int i = 0; i < num_of_elem; i++) {
		cout << "===========  a = " << a[i] << ", b = " << b[i] << " ===========" << endl;
		serial_time(a[i], b[i], eps[i], initial_step);
		for (int num_of_threads = 2; num_of_threads < 17; num_of_threads + 2) {
			parallel_time(a[i], b[i], eps[i], initial_step, num_of_threads);
		}
		cout << endl;
	}
}