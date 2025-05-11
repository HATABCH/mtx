#ifndef MATRIX_OPERATIONS_H
#define MATRIX_OPERATIONS_H

#include "matrix.h"

#define EPSILON 1e-10
#define MAX_ITERATIONS 1000
#define MTX_MAX_EXP_ITERATIONS 20

int mtx_add(matrix_t *m1, const matrix_t *m2);
int mtx_sub(matrix_t *m1, const matrix_t *m2);
int mtx_mul(matrix_t *m1, const matrix_t *m2);
int mtx_mul3(matrix_t *result, const matrix_t *m1, const matrix_t *m2);
int mtx_div(matrix_t *m1, const matrix_t *m2);

int mtx_add_scaled(matrix_t *m1, const matrix_t *m2, double scale);
int mtx_scale(matrix_t *m, double scale);

double mtx_norm(const matrix_t *m);

int mtx_inverse(const matrix_t *m, matrix_t *inv);
int mtx_gauss_elimination(matrix_t *m);
int mtx_exp(const matrix_t *m, matrix_t *result);

#endif // MATRIX_OPERATIONS_H
