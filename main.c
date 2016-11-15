#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdbool.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct Tableau {
	uint32_t M;
	uint32_t N;
	double *a;
	double *b;
} Tableau_t;

inline double f(const double x[4]) {
	return 2.4 * x[0] + 2.7 * x[1] + 13.8 * x[2] + 2.75 * x[3];
}

int32_t triangle(double *mat, double *b, uint32_t N, uint32_t M);

int32_t print(const double *mat, uint32_t N, uint32_t M) {
	for (uint32_t i = 0; i < M; ++i) {
		for (uint32_t j = 0; j < N; ++j) {
			printf("%.2f\t", mat[i * N + j]);
		}
		printf("\n");
	}

	printf("\n");

	return 0;
}

double * simplex_max(Tableau_t a);

int32_t main() {
	Tableau_t a;
	const uint32_t dim = 4;
	a.M = dim;
	a.N = dim;
	//a.a = (double *)malloc(a.N * a.M * sizeof(double));
	double mat[4 * 4] = { 1.01, 1.01, 9.45, 0.95,
                          0.25, 1.0 / 6.0, 0, 0,
	                      0.0, 0.0, 1.0 / 30.0, 4.0,
	                      0.0, 0.0, 0.0, 1.0 };
	double b[4] = { 140.0, 21.0, 16.0, 15.0 };
	a.a = mat;
	a.b = b;
	double *x = simplex_max(a);

	system("pause");

	return 0;
}

int32_t triangle(double *mat, double *b, const uint32_t N, const uint32_t M) {
	for (uint32_t i = 0; i < M; ++i) {
		double *denumerator = mat + i * N + i;
		for (uint32_t k = 0; k < M; ++k) {
			if (k == i) {
				const double tmp = 1.0 / *denumerator;
				for (uint32_t j = i; j < N; ++j) {
					mat[i * N + j] *= tmp;
				}
				b[i] *= tmp;
			} else {
				const double tmp = mat[k * N + i] / *denumerator;
				for (uint32_t j = i; j < N; ++j) {
					mat[k * N + j] -= tmp * mat[i * N + j];
				}
			}

		}
	}

	return 0;
}

double * simplex_max(Tableau_t a) {
	const double T = 1.0e+10;
	const uint32_t N = a.N + a.M;
	const uint32_t M = a.M;
	double *a_ext = (double *)calloc(N * M, sizeof(double));
	double *b = (double *)malloc(M * sizeof(double));
	memcpy(b, a.b, M * sizeof(double));
	for (uint32_t i = 0; i < M; ++i) {
		memcpy(a_ext + i * N, a.a + i * a.N, a.N * sizeof(double));
		*(a_ext + i * N + M + i) = 1.0;
	}
	//print(a_ext, N, M);
	triangle(a_ext, b, N, M);
	//print(a_ext, N, M);

	return 0;
}