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

inline double f(const double *c, const double *x, const uint32_t N);
int32_t triangle(double *mat, double *b, const uint32_t N, const uint32_t M);
void print_ab(const double *mat, const double *b, const uint32_t M, const uint32_t N);
void print_x(const double *x, const uint32_t N);
void replaceColMat(double *mat, const uint32_t i, const uint32_t j, const uint32_t M, const uint32_t N);
void reorderX(double *x, const uint32_t *reoder, const uint32_t M);
int32_t gauss(const double *a, const double *b, double *x, const uint32_t start, const uint32_t end, const uint32_t N);
int32_t findX(const double *a, const double *b, double *x, const uint32_t M, const uint32_t K);
int32_t simplex_max(const Tableau_t *table, double *x);

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

	double *x = (double *)calloc(a.M, sizeof(double));
	simplex_max(&a, x);

	system("pause");

	return 0;
}

inline double f(const double *c, const double *x, const uint32_t N) {
	double res = 0.0;
	for (uint32_t i = 0; i < N; ++i) {
		res += c[i] * x[i];
	}
	return res;
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
			}
			else {
				const double tmp = mat[k * N + i] / *denumerator;
				for (uint32_t j = i; j < N; ++j) {
					mat[k * N + j] -= tmp * mat[i * N + j];
				}
				b[k] -= tmp * b[i];
			}

		}
	}

	return 0;
}

void print_ab(const double *mat, const double *b, const uint32_t M, const uint32_t N) {
	for (uint32_t i = 0; i < M; ++i) {
		for (uint32_t j = 0; j < N; ++j) {
			printf("%.2f\t", mat[i * N + j]);
		}
		printf(" | %.2f\n", b[i]);
	}
	printf("\n");

	return;
}

void print_x(const double *x, const uint32_t N) {
	for (uint32_t i = 0; i < N; ++i) {
		printf("x[%i] = %.2f\n", i, x[i]);
	}
	printf("\n");

	return;
}

void replaceColMat(double *mat, const uint32_t i, const uint32_t j, const uint32_t M, const uint32_t N) {
	for (uint32_t k = 0; k < M; ++k) {
		double *mat_ki = mat + k * N + i;
		double *mat_kj = mat + k * N + j;
		const double temp = *mat_ki;
		*(mat_ki) = *(mat_kj);
		*(mat_kj) = temp;
	}

	return;
}

void reorderX(double *x, const uint32_t *reoder, const uint32_t M) {
	for (uint32_t i = M; i > 0; --i) {
		const uint32_t j = reoder[i - 1];
		double *x_i = x + i - 1;
		double *x_j = x + j;

		const double temp = *x_i;
		*x_i = *x_j;
		*x_j = temp;
	}

	return;
}

int32_t gauss(const double *a, const double *b, double *x, const uint32_t start, const uint32_t end, const uint32_t N) {
	const uint32_t M = end < N ? N - end : N - start;
	double *A = (double *)malloc(M * M * sizeof(double));
	for (uint32_t i = 0; i < M; ++i) {
		memcpy(A + i * M, a + i * N + start, M * sizeof(double));
	}
	double *rhs = (double *)malloc(M * sizeof(double));
	memcpy(rhs, b, M * sizeof(double));
	double *x_ = x + start;

	double max;
	double factor;
	uint32_t* replaceCol = (uint32_t *)malloc(M * sizeof(uint32_t));
	for (uint32_t i = 0; i < M; ++i) {
		replaceCol[i] = i;
		max = fabs(*(A + i * M + i));
		for (uint32_t k = i + 1; k < M; ++k) {
			if (fabs(*(A + i * M + k)) > max) {
				max = fabs(*(A + i * M + k));
				replaceCol[i] = k;
			}
		}

		if (replaceCol[i] != i) {
			replaceColMat(A, i, replaceCol[i], M, M);
		}

		for (uint32_t k = i + 1; k < M; ++k) {
			double *a_ki = A + k * M + i;
			if (*a_ki == 0.0)
				continue;

			factor = *a_ki / *(A + i * M + i);
			*a_ki = 0.0;
			for (uint32_t j = i + 1; j < M; ++j)
				*(A + k * M + j) -= *(A + i * M + j) * factor;
			rhs[k] -= rhs[i] * factor;
		}
	}

	int32_t k;
	for (int32_t i = M - 1; i > 0; --i) {
		x_[i] = rhs[i] / A[i * M + i];
		k = 1;
		while (i - k >= 0) {
			rhs[i - k] -= A[(i - k) * M + i] * x_[i];
			k++;
		}
	}
	x_[0] = rhs[0] / A[0];

	reorderX(x_, replaceCol, M);

	free(A);
	free(rhs);
	free(replaceCol);
	return 0;
}

int32_t simplex_max(const Tableau_t *table, double *x) {
	const double T = -1.0e+10;
	const uint32_t N = table->N;
	const uint32_t M = table->M;
	const uint32_t K = M + N;
	double *a = (double *)calloc(M * K, sizeof(double));
	double *b = (double *)malloc(M * sizeof(double));
	double *c = (double *)malloc(K * sizeof(double));
	double *delta = (double *)calloc(M, sizeof(double));
	double *x_ = (double *)calloc(K, sizeof(double));
	//uint32_t *col_index = (uint32_t *)malloc(K * sizeof(uint32_t));
	memcpy(x_, x, M * sizeof(double));
	memcpy(b, table->b, M * sizeof(double));
	memcpy(c, table->c, M * sizeof(double));
	for (uint32_t i = 0; i < M; ++i) {
		memcpy(a + i * K, table->a + i * N, N * sizeof(double));
		a[i * K + i + M] = 1.0;
		//col_index[i] = i;
	}
	for (uint32_t i = M; i < K; ++i) {
		c[i] = T;
		//col_index[i] = i;
	}
	//print_ab(a, b, M, K);

	gauss(a, b, x_, M, K, K);

	for (uint32_t z = 0; z < 20; ++z) {
		for (uint32_t i = 0; i < M; ++i) {
			double z = 0.0;
			const uint32_t l = M + i;
			for (uint32_t k = 0; k < M; ++k) {
				z += c[k + M] * a[k * K + l];
			}
			delta[i] = c[i] - z;
		}
		/*
		for (uint32_t i = 0; i < M; ++i) {
			printf("delta[%i] = %.2f\n", i, delta[i]);
		}
		printf("\n");
		*/
		if (delta[0] <= 0 && delta[1] <= 0 && delta[2] <= 0 && delta[3] <= 0) break;

		double max_delta = 0.0;
		uint32_t col_r = 0;
		for (uint32_t i = 0; i < M; ++i) {
			const double tmp = delta[i];
			if (tmp > max_delta) {
				max_delta = tmp;
				col_r = i;
			}
		}

		if (max_delta != 0.0) {
			uint32_t col_s = 0;
			uint32_t row_s = 0;
			const double max = 1.0e+10;
			double min_delta = max;
			for (uint32_t i = 0; i < M; ++i) {
				if (a[i * K + col_r] <= 0) continue;
				const double tmp = x_[i + M] / a[i * K + col_r];
				if (tmp <= min_delta) {
					min_delta = tmp;
					row_s = i;
					col_s = i + M;
				}
			}

			if (min_delta != max) {
				const double tmp = 1.0 / a[row_s * K + col_r];
				x_[col_s] *= tmp;
				for (uint32_t i = 0; i < K; ++i) {
					a[row_s * K + i] *= tmp;
				}

				for (uint32_t i = 0; i < M; ++i) {
					if (i == row_s) continue;
					const double tmp2 = a[i * K + col_r] / a[row_s * K + col_r];
					for (uint32_t j = 0; j < K; ++j) {
						a[i * K + j] -= a[row_s * K + j] * tmp2;
					}
					x_[i + M] -= x_[col_s] * tmp2;
				}

				replaceColMat(a, col_r, col_s, M, K);
				//replaceColMat(col_index, col_r, col_s, 1, K);
				replaceColMat(c, col_r, col_s, 1, K);
			} else {
				break;
			}
		} else {
			printf("No Positive!\n");
			break;
		}
	}
	print_x(x_, K);
	//printf("F(x) = %f\n", f(c, x_, K));

	free(a);
	free(b);
	free(c);
	free(x_);
	free(delta);

	//memcpy(x, x_, M * sizeof(double));

	return 0;
}