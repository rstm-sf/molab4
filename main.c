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
void print_ab(const double *mat, const double *b, const uint32_t M, const uint32_t N);
void print_x(const double *x, const uint32_t N);
void replaceColMat(double *mat, const uint32_t i, const uint32_t j, const uint32_t M, const uint32_t N);
void reorderX(double *x, const uint32_t *reoder, const uint32_t M);
int32_t simplex_max(const Tableau_t *table, double *x);

int32_t main() {
	Tableau_t table;
	const uint32_t dim = 4;
	table.M = dim;
	table.N = dim;
	double mat[4 * 4] = { 1.01, 1.01, 9.45, 0.95,
	                      0.25, 1.0 / 6.0, 0, 0,
	                      0.0, 0.0, 1.0 / 30.0, 4.0,
	                      0.0, 0.0, 0.0, 1.0 };
	double b[4] = { 140.0, 21.0, 16.0, 15.0 };
	double c[4] = { 2.4, 2.7, 13.8, 2.75 };
	table.a = mat;
	table.b = b;
	table.c = c;

	printf("F(x) = ");
	for (uint32_t i = 0; i < table.M; ++i) {
		printf("+ (%.2f) * x[%" PRIu32 "] ", table.c[i], i);
	}
	printf("\n");
	print_ab(table.a, table.b, table.M, table.N);

	double *x = (double *)calloc(table.M, sizeof(double));
	simplex_max(&table, x);

	print_x(x, table.M);
	printf("F(x) = %f\n", f(table.c, x, table.M));

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

void print_ab(const double *mat, const double *b, const uint32_t M, const uint32_t N) {
	for (uint32_t i = 0; i < M; ++i) {
		for (uint32_t j = 0; j < N; ++j) {
			printf("%.2f\t", mat[i * N + j]);
		}
		printf("<=\t%.2f\n", b[i]);
	}
	printf("\n");

	return;
}

void print_x(const double *x, const uint32_t N) {
	for (uint32_t i = 0; i < N; ++i) {
		printf("x[%" PRIu32 "] = %f\n", i, x[i]);
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
	double *tmp = (double *)malloc(M * sizeof(double));
	memcpy(tmp, x, M * sizeof(double));

	for (uint32_t i = 0; i < M; ++i) {
		*(x) = *(tmp + *(reoder + i));
		x++;
	}

	free(tmp);

	return;
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
	uint32_t *col_index = (uint32_t *)malloc(K * sizeof(uint32_t));

	memcpy(x_, x, M * sizeof(double));
	memcpy(b, table->b, M * sizeof(double));
	memcpy(c, table->c, M * sizeof(double));
	for (uint32_t i = 0; i < M; ++i) {
		memcpy(a + i * K, table->a + i * N, N * sizeof(double));
		a[i * K + i + N] = 1.0;
		col_index[i] = i;
	}
	for (uint32_t i = M; i < K; ++i) {
		if (i >= N) {
			c[i] = T;
		}
		col_index[i] = i;
	}

	for (uint32_t i = K - M; i < K; ++i) {
		x_[i] = b[i - M];
	}

	uint32_t z;
	const uint32_t max_iter = 20;
	for (z = 0; z < max_iter; ++z) {
		bool isPositiveDelta = false;

		for (uint32_t i = 0; i < M; ++i) {
			double z = 0.0;
			const uint32_t l = M + i;
			for (uint32_t k = 0; k < M; ++k) {
				z += c[k + M] * a[k * K + l];
			}
			delta[i] = c[i] - z;

			if (delta[i] > 0.0) {
				isPositiveDelta = true;
			}
		}

		if (!isPositiveDelta) break;

		double max_delta = 0.0;
		uint32_t col_r = 0;
		for (uint32_t i = 0; i < M; ++i) {
			const double tmp = delta[i];
			if (tmp > max_delta) {
				max_delta = tmp;
				col_r = i;
			}
		}

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

		const double tmp = 1.0 / a[row_s * K + col_r];
		x_[col_s] *= tmp;
		for (uint32_t i = 0; i < K; ++i) {
			a[row_s * K + i] *= tmp;
		}

		for (uint32_t i = 0; i < M; ++i) {
			if (i == row_s) continue;
			const double factor = a[i * K + col_r] / a[row_s * K + col_r];
			for (uint32_t j = 0; j < K; ++j) {
				a[i * K + j] -= a[row_s * K + j] * factor;
			}
			x_[i + M] -= x_[col_s] * factor;
		}

		replaceColMat(a, col_r, col_s, M, K);
		replaceColMat(c, col_r, col_s, 1, K);
		const uint32_t tmp2 = col_index[col_r];
		col_index[col_r] = col_index[col_s];
		col_index[col_s] = tmp2;
	}
	if (z == max_iter) {
		printf("Warning: MAX ITER!\n");
	}
	printf("ITER = %" PRId32 "\n\n", z);

	reorderX(x_, col_index, K);
	memcpy(x, x_, M * sizeof(double));

	free(a);
	free(b);
	free(c);
	free(x_);
	free(delta);

	return 0;
}