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
	double *c;
} Tableau_t;

inline double f(const double *c, const double *x, const uint32_t N) {
	double res = 0.0;
	for (uint32_t i = 0; i < N; ++i) {
		res += c[i] * x[i];
	}
	return res;
}

int32_t triangle(double *mat, double *b, uint32_t N, uint32_t M);
int32_t gauss(double *a, double *b, double *x, uint32_t M, uint32_t N);

int32_t print(const double *mat, const double *b, uint32_t N, uint32_t M) {
	for (uint32_t i = 0; i < M; ++i) {
		for (uint32_t j = 0; j < N; ++j) {
			printf("%.2f\t", mat[i * N + j]);
		}
		printf("%.2f\n", b[i]);
	}
	printf("\n");

	return 0;
}

double * simplex_max(Tableau_t *a);

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
	double c[4] = { 2.4, 2.7, 13.8, 2.75 };
	a.a = mat;
	a.b = b;
	a.c = c;
	double *x = simplex_max(&a);

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

int32_t gauss(const double *a, const double *b, double *x, const uint32_t M, const uint32_t N) {
	double *A = (double *)malloc(M * M * sizeof(double));
	for (uint32_t i = 0; i < M; ++i) {
		for (uint32_t j = 0; j < M; ++j) {
			A[i * M + j] = a[i * N + j + M];
		}
	}
	double *rhs = (double *)malloc(M * sizeof(double));
	memcpy(rhs, b, M * sizeof(double));

	double temp, max;
	double factor;
	uint32_t* replaceCol = (uint32_t *)malloc(M * sizeof(uint32_t));
	for (uint32_t i = 0; i < M; ++i) {
		replaceCol[i] = i;
		max = fabs(A[replaceCol[i] + i * M]);
		for (uint32_t k = i + 1; k < M; ++k) {
			if (fabs(A[i * M + k]) > max) {
				max = fabs(A[i * M + k]);
				replaceCol[i] = k;
			}
		}
		if (replaceCol[i] != i) {
			for (uint32_t k = 0; k < M; ++k) {
				temp = A[k * M + i];
				A[k * M + i] = A[k * M + replaceCol[i]];
				A[k * M + replaceCol[i]] = temp;
			}
		}
		for (uint32_t k = i + 1; k < M; ++k) {
			if (A[k * M + i] == 0)
				continue;

			factor = A[k * M + i] / A[i * M + i];
			A[k * M + i] = 0;
			for (uint32_t j = i + 1; j < M; ++j)
				A[k * M + j] = A[k * M + j] - A[i * M + j] * factor;
			rhs[k] = rhs[k] - rhs[i] * factor;
		}
	}

	uint32_t k;
	for (uint32_t i = M - 1; i > 0; --i) {
		x[i] = rhs[i] / A[i * M + i];
		k = 1;
		while (i - k >= 0) {
			rhs[i - k] = rhs[i - k] - A[(i - k) * M + i] * x[i];
			k++;
		}
	}
	x[0] = rhs[0] / A[0];

	for (uint32_t i = M - 1; i >= 0; --i) {
		temp = x[i];
		x[i] = x[replaceCol[i]];
		x[replaceCol[i]] = temp;
	}

	free(A);
	free(rhs);
	free(replaceCol);
	return 0;
}

double * simplex_max(const Tableau_t *a) {
	const double T = -1.0e+10;
	const uint32_t N = a->N + a->M;
	const uint32_t M = a->M;
	double *a_ext = (double *)calloc(N * M, sizeof(double));
	double *b = (double *)malloc(M * sizeof(double));
	double *c = (double *)malloc(N * sizeof(double));
	memcpy(b, a->b, M * sizeof(double));
	memcpy(c, a->c, M * sizeof(double));
	for (uint32_t i = M; i < N; ++i) {
		c[i] = T;
	}
	for (uint32_t i = 0; i < M; ++i) {
		memcpy(a_ext + i * N, a->a + i * a->N, a->N * sizeof(double));
		*(a_ext + i * N + M + i) = 1.0;
	}
	//print(a_ext, b, N, M);
	triangle(a_ext, b, N, M);
	print(a_ext, b, N, M);
	double *delta = (double *)calloc(M, sizeof(double));
	double *x = (double *)calloc(N, sizeof(double));

	for (uint32_t i = 0; i < M; ++i) {
		double z = 0.0;
		for (uint32_t k = 0; k < M; ++k) {
			z += c[k] * a_ext[k * N + i];
		}
		delta[i] = c[i + M] - z;
	}
	/*
	for (uint32_t i = 0; i < M; ++i) {
		printf("%.2f\t", delta[i]);
	}
	printf("\n");
	*/
	free(a_ext);
	free(b);
	free(c);
	free(delta);
	free(x);

	return 0;
}