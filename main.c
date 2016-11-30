#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdbool.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <assert.h>

typedef enum Sign {
	NOT_GREATER_THAN = 1,
	EQUALITY = 0,
	NOT_LESS_THAN = -1
} Sign_t;

typedef enum Extreme {
	MIN = -1,
	MAX = 1
} Extreme_t;

typedef struct Tableau {
	uint32_t M;
	uint32_t N;
	double *a;
	double *b;
	double *c;
	Sign_t *sign;
	Extreme_t extr;
} Tableau_t;

inline double f(const double *c, const double *x, const uint32_t N);
void print_a_task(const double *mat, const uint32_t M, const uint32_t N);
void print_x(const double *x, const uint32_t N);
void copy_c_cond_extr(double *dst_c, const double *src_c, const uint32_t N, const Extreme_t extreme);
int32_t simplex(const Tableau_t *table, double *x);

int32_t main() {
	double b[4] = { 140.0, 21.0, 16.0, 15.0 };
	double c[4] = { 2.4, 2.7, 13.8, 2.75 };
	Extreme_t extr = MIN;

	Tableau_t table;
	table.M = table.N = 4;
	table.extr = extr;
	if (extr == MAX) {
		double mat[4 * 4] = { 1.01, 1.01,      9.45,       0.95,
	                          0.2,  1.0 / 6.0, 0.0,        0.0,
	                          0.0,  0.0,       1.0 / 30.0, 4.0,
	                          0.0,  0.0,       0.0,        1.0 };
		Sign_t sign[4] = { NOT_GREATER_THAN, NOT_GREATER_THAN, NOT_GREATER_THAN, NOT_GREATER_THAN };

		table.a = mat;
		table.b = b;
		table.c = c;
		table.sign = sign;
	} else {
		double mat[4 * 4] = { 1.01, 0.2,       0.0,        0.0,
	                          1.01, 1.0 / 6.0, 0.0,        0.0,
	                          9.45, 0.0,       1.0 / 30.0, 0.0,
	                          0.95, 0.0,       4.0,        1.0 };
		Sign_t sign[4] = { NOT_LESS_THAN, NOT_LESS_THAN, NOT_LESS_THAN, NOT_LESS_THAN };

		table.a = mat;
		table.b = c;
		table.c = b;
		table.sign = sign;
	}

	double *x = (double *)calloc(table.N, sizeof(double));
	simplex(&table, x);

	print_x(x, table.M);
	printf("F(x) = %f\n", f(table.c, x, table.N));

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

void print_a_task(const double *mat, const uint32_t M, const uint32_t N) {
	const double *mat_ = mat;
	for (uint32_t i = 0; i < M; ++i) {
		for (uint32_t j = 0; j < N; ++j) {
			printf("%.4f\t", *mat_++);
		}
		printf("\n");
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

// if extreme = MIN then c_i = -ci
void copy_c_cond_extr(double *dst_c, const double *src_c, const uint32_t N, const Extreme_t extreme) {
	double       *p = dst_c;
	const double *q = src_c;
	for (uint32_t i = 0; i < N; ++i) {
		*p++ = *q++ * extreme;
	}
	return;
}

int32_t simplex(const Tableau_t *table, double *x) {
	static const double sufficiently_large_number = 1.0e+10;
	const uint32_t N = table->N;
	const uint32_t M = table->M;
	const uint32_t K = M; // in the number inequality
	const uint32_t N_ext = N + K + M;
	const uint32_t offset_b = N_ext - M;
	const double *b = table->b;

	for (uint32_t i = 0; i < N; ++i) {
		assert(b[i] >= 0);
	}

	double *a_ = (double *)calloc(M * N_ext, sizeof(double));
	double *c_ = (double *)malloc(N_ext * sizeof(double));
	double *x_ = (double *)calloc(N_ext, sizeof(double));
	// col_index - the room basic and not basic columns
	// col_index[0 .. N_ext - M - 1] - not basic; col_index[N_ext - M .. N_ext - 1] - basic
	uint32_t *col_index = (uint32_t *)malloc(N_ext * sizeof(uint32_t));
	uint32_t *col_index_nb = col_index;
	uint32_t *col_index_b = col_index + offset_b;

	memcpy(x_, x, N * sizeof(double));
	copy_c_cond_extr(c_, table->c, N, table->extr);
	for (uint32_t i = 0; i < M; ++i) {
		memcpy(a_ + i * N_ext, table->a + i * N, N * sizeof(double));
		a_[i * N_ext + i + N] = (double)(table->sign[i]);
		a_[i * N_ext + i + offset_b] = 1.0;

		col_index[i] = i;
	}
	for (uint32_t i = M; i < N_ext; ++i) {
		if (i >= N && i < offset_b) {
			c_[i] = 0.0;
		} else {
			c_[i] = -sufficiently_large_number;
		}

		col_index[i] = i;
	}

	for (uint32_t i = offset_b; i < N_ext; ++i) {
		x_[i] = b[i - offset_b];
	}

	uint32_t iter;
	const uint32_t max_iter = N_ext + 1;
	for (iter = 0; iter < max_iter; ++iter) {
		printf("ITER = %" PRId32 "\n", iter);
		print_a_task(a_, M, N_ext);

		uint32_t index_r = 0;

		bool isExitDelta = true;
		double max_delta = 0.0;
		for (uint32_t i = 0; i < offset_b; ++i) {
			double z = 0.0;
			const uint32_t ii = col_index_nb[i];
			for (uint32_t k = 0; k < M; ++k) {
				const uint32_t kk = col_index_b[k];
				z += c_[kk] * a_[k * N_ext + ii];
			}
			const double delta_i = c_[ii] - z;

			if (delta_i > max_delta) {
				isExitDelta = false;
				max_delta = delta_i;
				index_r = i;
			}
		}

		if (isExitDelta) {
			break;
		}

		const uint32_t col_r = col_index_nb[index_r];

		uint32_t row_s = 0;

		double min_delta = DBL_MAX;
		for (uint32_t i = 0; i < M; ++i) {
			const double a_ir = a_[i * N_ext + col_r];
			if (a_ir <= 0) continue;
			const double tmp = x_[col_index_b[i]] / a_ir;
			if (tmp < min_delta) {
				min_delta = tmp;
				row_s = i;
			} else if (tmp == min_delta) {
				row_s = row_s > i ? row_s : i;
			}
		}

		const uint32_t col_s = col_index_b[row_s];
		const double x_r = min_delta;
		double *a_s = a_ + row_s * N_ext;

		const double tmp = 1.0 / *(a_s + col_r);
		for (uint32_t i = 0; i < N_ext; ++i) {
			*(a_s + i) *= tmp;
		}

		for (uint32_t i = 0; i < M; ++i) {
			if (i == row_s) continue;
			double *a_i = a_ + i * N_ext;
			const double factor = *(a_i + col_r); // a_[s,r] = 1
			x_[col_index_b[i]] -= factor * x_r;
			for (uint32_t j = 0; j < N_ext; ++j) {
				*(a_i + j) -= *(a_s + j) * factor;
			}
		}

		col_index_b[row_s] = col_r;
		col_index_nb[index_r] = col_s;
		x_[col_s] = 0.0;
		x_[col_r] = x_r;
	}
	if (iter == max_iter) {
		printf("\nWarning: MAX ITER!\n");
	}
	printf("ITERs = %" PRId32 "\n", iter);
	print_a_task(a_, M, N_ext);

	memcpy(x, x_, N * sizeof(double));

	free(a_);
	free(c_);
	free(x_);

	return 0;
}