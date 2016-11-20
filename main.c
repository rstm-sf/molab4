#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdbool.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

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
int32_t simplex_max(const Tableau_t *table, double *x);

int32_t main() {
	double mat[4 * 4] = { 1.01, 1.01,      9.45,       0.95,
	                      0.25, 1.0 / 6.0, 0.0,        0.0,
	                      0.0,  0.0,       1.0 / 30.0, 4.0,
	                      0.0,  0.0,       0.0,        1.0 };
	double b[4] = { 140.0, 21.0, 16.0, 15.0 };
	double c[4] = { 2.4, 2.7, 13.8, 2.75 };

	Tableau_t table;
	table.M = 4;
	table.N = 4;
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

	free(x);
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

int32_t simplex_max(const Tableau_t *table, double *x) {
	const double T = 1.0e+10;
	const uint32_t N = table->N;
	const uint32_t M = table->M;
	const uint32_t K = M + N;
	const double *b = table->b;

	double *a = (double *)calloc(M * K, sizeof(double));
	double *c = (double *)malloc(K * sizeof(double));
	double *x_ = (double *)calloc(K, sizeof(double));
	double *delta = (double *)calloc(K - M, sizeof(double));
	// col_index - the room basic and not basic columns
	// col_index[0 .. K - M] - not basic; col_index[K - M .. K - 1] - basic
	uint32_t *col_index = (uint32_t *)malloc(K * sizeof(uint32_t));

	memcpy(x_, x, M * sizeof(double));
	memcpy(c, table->c, M * sizeof(double));
	for (uint32_t i = 0; i < M; ++i) {
		memcpy(a + i * K, table->a + i * N, N * sizeof(double));
		a[i * K + i + N] = 1.0;

		col_index[i] = i;
	}
	for (uint32_t i = M; i < K; ++i) {
		if (i >= N) {
			c[i] = -T;
		}
		col_index[i] = i;
	}

	for (uint32_t i = K - M; i < K; ++i) {
		x_[i] = b[i - K + M];
	}

	uint32_t iter;
	const uint32_t max_iter = 20;
	for (iter = 0; iter < max_iter; ++iter) {
		bool isPositiveDelta = false;
		//print_ab(a, x_ + M, M, K);

		for (uint32_t i = 0; i < M; ++i) {
			double z = 0.0;
			const uint32_t ii = col_index[i];
			for (uint32_t k = 0; k < M; ++k) {
				const uint32_t kk = col_index[k + K - M];
				z += c[kk] * a[k * K + ii];
			}
			delta[i] = c[ii] - z;

			if (delta[i] > 0.0) {
				isPositiveDelta = true;
			}
		}

		if (!isPositiveDelta) break;

		double max_delta = 0.0;
		uint32_t index_r = 0;
		for (uint32_t i = 0; i < M; ++i) {
			const double tmp = delta[i];
			if (tmp > max_delta) {
				max_delta = tmp;
				index_r = i;
			}
		}
		const uint32_t col_r = col_index[index_r];

		uint32_t row_s = 0;
		double min_delta = DBL_MAX;
		for (uint32_t i = 0; i < M; ++i) {
			const double a_ir = a[i * K + col_r];
			if (a_ir <= 0) continue;
			const double tmp = x_[col_index[i + K - M]] / a_ir;
			if (tmp <= min_delta) {
				min_delta = tmp;
				row_s = i;
			}
		}
		const uint32_t col_s = col_index[row_s + M];
		double *a_s = a + row_s * K;

		const double tmp = 1.0 / *(a + row_s * K + col_r);
		x_[col_r] = x_[col_s] * tmp;
		x_[col_s] = 0.0;
		const double x_r = x_[col_r];
		for (uint32_t i = 0; i < K; ++i) {
			*(a_s + i) *= tmp;
		}

		for (uint32_t i = 0; i < M; ++i) {
			if (i == row_s) continue;
			double *a_i = a + i * K;
			const double factor = *(a_i + col_r); // a[s,r] = 1
			x_[col_index[i + M]] -= factor * x_r;
			for (uint32_t j = 0; j < K; ++j) {
				*(a_i + j) -= *(a_s + j) * factor;
			}
		}

		col_index[row_s + M] = col_r;
		col_index[index_r] = col_s;
	}
	if (iter == max_iter) {
		printf("Warning: MAX ITER!\n");
	}
	printf("ITER = %" PRId32 "\n\n", iter);

	memcpy(x, x_, M * sizeof(double));

	free(a);
	free(c);
	free(x_);
	free(delta);

	return 0;
}