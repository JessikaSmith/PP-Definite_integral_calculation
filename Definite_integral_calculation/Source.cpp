#include <iostream>
#include <omp.h>
#include <stdlib.h>  

using namespace std;

static long npoints = 100;
double eps = 0.001;

int main() {
	float a, b;

	for (int i = 0; i <= npoints; i++) {
		double J = 1 / 4*(2*(b - a) / (a*b) + sin(2 / b) - sin(2 / a));
	}
	cin >> a;
	cout << 'a';
}