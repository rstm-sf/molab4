#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdbool.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <assert.h>

typedef struct Tableau {
	uint32_t M;
	uint32_t N;
	double *a;
	double *b;
	double *c;
	void(*print)(const struct Tableau *);
} Tableau_t;

inline void print_task(const Tableau_t *table);
inline double f(const double *c, const double *x, const uint32_t N);
void print_ab(const double *mat, const double *b, const uint32_t M, const uint32_t N);
void print_x(const double *x, const uint32_t N);
int32_t simplex_max(const Tableau_t *table, double *x);
int32_t simplex_min(const Tableau_t *table, double *x);

int32_t main() {
	/*
	double mat[4 * 4] = { 1.01, 1.01,      9.45,       0.95,
	                      0.25, 1.0 / 6.0, 0.0,        0.0,
	                      0.0,  0.0,       1.0 / 30.0, 4.0,
	                      0.0,  0.0,       0.0,        1.0 };
	double b[4] = { 140.0, 21.0, 16.0, 15.0 };
	double c[4] = { 2.4, 2.7, 13.8, 2.75 };
	*/
	double mat[4 * 4] = { 1.01, 0.25,      0.0,        0.0,
	                      1.01, 1.0 / 6.0, 0.0,        0.0,
	                      9.45, 0.0,       1.0 / 30.0, 0.0,
	                      0.95, 0.0,       4.0,        1.0 };
	double b[4] = { 2.4, 2.7, 13.8, 2.75 };
	double c[4] = { 140.0, 21.0, 16.0, 15.0 };

	Tableau_t table;
	table.M = 4;
	table.N = 4;
	table.a = mat;
	table.b = b;
	table.c = c;
	table.print = print_task;
	table.print(&table);

	double *x = (double *)calloc(table.M, sizeof(double));
	simplex_min(&table, x);

	print_x(x, table.M);
	printf("F(x) = %f\n", f(table.c, x, table.M));

	system("pause");

	free(x);
	return 0;
}

inline void print_task(const Tableau_t *table) {
	printf("F(x) = ");
	const uint32_t M = table->M;
	const double *c = table->c;
	for (uint32_t i = 0; i < M; ++i) {
		printf("+ (%.2f) * x[%" PRIu32 "] ", c[i], i);
	}
	printf("\n");
	print_ab(table->a, table->b, table->M, table->N);

	return;
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
	const uint32_t N_ext = M + N;
	const uint32_t offset_b = N_ext - M;
	const double *b = table->b;

	for (uint32_t i = 0; i < M; ++i) {
		assert(b[i] >= 0);
	}

	double *a = (double *)calloc(M * N_ext, sizeof(double));
	double *c = (double *)malloc(N_ext * sizeof(double));
	double *x_ = (double *)calloc(N_ext, sizeof(double));
	// col_index - the room basic and not basic columns
	// col_index[0 .. N_ext - M - 1] - not basic; col_index[N_ext - M .. N_ext - 1] - basic
	uint32_t *col_index = (uint32_t *)malloc(N_ext * sizeof(uint32_t));
	uint32_t *col_index_nb = col_index;
	uint32_t *col_index_b  = col_index + offset_b;

	memcpy(x_, x, N * sizeof(double));
	memcpy(c, table->c, N * sizeof(double));
	for (uint32_t i = 0; i < M; ++i) {
		memcpy(a + i * N_ext, table->a + i * N, N * sizeof(double));
		a[i * N_ext + i + N] = 1.0;

		col_index[i] = i;
	}
	for (uint32_t i = M; i < N_ext; ++i) {
		if (i >= N) {
			c[i] = -T;
		}
		col_index[i] = i;
	}

	for (uint32_t i = offset_b; i < N_ext; ++i) {
		x_[i] = b[i - offset_b];
	}

	uint32_t iter;
	const uint32_t max_iter = N_ext + 1;
	for (iter = 0; iter < max_iter; ++iter) {
		uint32_t index_r = 0;

		bool isPositiveDelta = false;
		double max_delta = 0.0;
		for (uint32_t i = 0; i < offset_b; ++i) {
			double z = 0.0;
			const uint32_t ii = col_index_nb[i];
			for (uint32_t k = 0; k < M; ++k) {
				const uint32_t kk = col_index_b[k];
				z += c[kk] * a[k * N_ext + ii];
			}
			const double delta_i = c[ii] - z;

			if (delta_i > max_delta) {
				isPositiveDelta = true;
				max_delta = delta_i;
				index_r = i;
			}
		}

		if (!isPositiveDelta) {
			break;
		}
		const uint32_t col_r = col_index_nb[index_r];

		uint32_t row_s = 0;
		double min_delta = DBL_MAX;
		for (uint32_t i = 0; i < M; ++i) {
			const double a_ir = a[i * N_ext + col_r];
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
		double *a_s = a + row_s * N_ext;

		const double tmp = 1.0 / *(a_s + col_r);
		x_[col_r] = x_[col_s] * tmp;
		const double x_r = x_[col_r];
		for (uint32_t i = 0; i < N_ext; ++i) {
			*(a_s + i) *= tmp;
		}

		for (uint32_t i = 0; i < M; ++i) {
			if (i == row_s) continue;
			double *a_i = a + i * N_ext;
			const double factor = *(a_i + col_r); // a[s,r] = 1
			x_[col_index_b[i]] -= factor * x_r;
			for (uint32_t j = 0; j < N_ext; ++j) {
				*(a_i + j) -= *(a_s + j) * factor;
			}
		}

		col_index_b[row_s] = col_r;
		col_index_nb[index_r] = col_s;
		x_[col_s] = 0.0;
	}
	if (iter == max_iter) {
		printf("Warning: MAX ITER!\n");
	}
	printf("ITER = %" PRId32 "\n\n", iter);

	memcpy(x, x_, M * sizeof(double));

	free(a);
	free(c);
	free(x_);

	return 0;
}

int32_t simplex_min(const Tableau_t *table, double *x) {
	const double T = 1.0e+10;
	const uint32_t N = table->N;
	const uint32_t M = table->M;
	const uint32_t N_ext = M + N + M;
	const uint32_t offset_b = N_ext - M;
	const double *b = table->b;

	for (uint32_t i = 0; i < M; ++i) {
		assert(b[i] >= 0);
	}

	double *a = (double *)calloc(M * N_ext, sizeof(double));
	double *c = (double *)malloc(N_ext * sizeof(double));
	double *x_ = (double *)calloc(N_ext, sizeof(double));
	// col_index - the room basic and not basic columns
	// col_index[0 .. N_ext - M - 1] - not basic; col_index[N_ext - M .. N_ext - 1] - basic
	uint32_t *col_index = (uint32_t *)malloc(N_ext * sizeof(uint32_t));
	uint32_t *col_index_nb = col_index;
	uint32_t *col_index_b = col_index + offset_b;

	memcpy(x_, x, N * sizeof(double));
	memcpy(c, table->c, N * sizeof(double));
	for (uint32_t i = 0; i < M; ++i) {
		memcpy(a + i * N_ext, table->a + i * N, N * sizeof(double));
		a[i * N_ext + i + N] = -1.0;
		a[i * N_ext + i + offset_b] = 1.0;

		col_index[i] = i;
	}
	for (uint32_t i = M; i < N_ext; ++i) {
		if (i >= N) {
			c[i] = T;
		}
		col_index[i] = i;
	}

	for (uint32_t i = offset_b; i < N_ext; ++i) {
		x_[i] = b[i - offset_b];
	}

	uint32_t iter;
	const uint32_t max_iter = N_ext + 1;
	for (iter = 0; iter < max_iter; ++iter) {
		uint32_t index_r = 0;

		bool isPositiveDelta = true;
		double max_delta = 0.0;
		for (uint32_t i = 0; i < offset_b; ++i) {
			double z = 0.0;
			const uint32_t ii = col_index_nb[i];
			for (uint32_t k = 0; k < M; ++k) {
				const uint32_t kk = col_index_b[k];
				z += c[kk] * a[k * N_ext + ii];
			}
			const double delta_i = c[ii] - z;

			if (delta_i < max_delta) {
				isPositiveDelta = false;
				max_delta = delta_i;
				index_r = i;
			}
		}

		if (isPositiveDelta) {
			break;
		}
		const uint32_t col_r = col_index_nb[index_r];

		uint32_t row_s = 0;
		double min_delta = DBL_MAX;
		for (uint32_t i = 0; i < M; ++i) {
			const double a_ir = a[i * N_ext + col_r];
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
		double *a_s = a + row_s * N_ext;

		const double tmp = 1.0 / *(a_s + col_r);
		x_[col_r] = x_[col_s] * tmp;
		const double x_r = x_[col_r];
		for (uint32_t i = 0; i < N_ext; ++i) {
			*(a_s + i) *= tmp;
		}

		for (uint32_t i = 0; i < M; ++i) {
			if (i == row_s) continue;
			double *a_i = a + i * N_ext;
			const double factor = *(a_i + col_r); // a[s,r] = 1
			x_[col_index_b[i]] -= factor * x_r;
			for (uint32_t j = 0; j < N_ext; ++j) {
				*(a_i + j) -= *(a_s + j) * factor;
			}
		}

		col_index_b[row_s] = col_r;
		col_index_nb[index_r] = col_s;
		x_[col_s] = 0.0;
	}
	if (iter == max_iter) {
		printf("Warning: MAX ITER!\n");
	}
	printf("ITER = %" PRId32 "\n\n", iter);

	memcpy(x, x_, M * sizeof(double));

	free(a);
	free(c);
	free(x_);

	return 0;
}