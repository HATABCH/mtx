#ifndef MATRIX_LIB_H
#define MATRIX_LIB_H

#include <stdbool.h>
#include <stddef.h>

struct matrix;
typedef struct matrix matrix_t;

#define MTX_EPSILON 1e-9
#define MTX_GAUSS_EPSILON 1e-12
#define MTX_MAX_EXP_ITERATIONS 100

matrix_t *mtx_alloc(size_t h, size_t w);
matrix_t *mtx_copy(const matrix_t *m);
void mtx_free(matrix_t *m);

double *mtx_ptr(matrix_t *m, size_t r, size_t c);
const double *mtx_cptr(const matrix_t *m, size_t r, size_t c);

size_t mtx_get_rows(const matrix_t *m);
size_t mtx_get_cols(const matrix_t *m);

void mtx_set_zero(matrix_t *m);
void mtx_set_id(matrix_t *m);

matrix_t *mtx_alloc_zero(size_t h, size_t w);
matrix_t *mtx_alloc_id(size_t h, size_t w);

int mtx_assign(matrix_t *m1_dst, const matrix_t *m2_src);

void mtx_move_assign(matrix_t *target, matrix_t *source_to_consume);

int mtx_read(matrix_t *m);
void mtx_print(const matrix_t *m);
void mtx_print_titled(const char *title, const matrix_t *m);

int mtx_add(matrix_t *m1_res, const matrix_t *m2_add);
int mtx_sub(matrix_t *m1_res, const matrix_t *m2_sub);

void mtx_smul(matrix_t *m, double scalar);
int mtx_sdiv(matrix_t *m, double scalar);

int mtx_add2(matrix_t *m_res, const matrix_t *m1, const matrix_t *m2);
int mtx_sub2(matrix_t *m_res, const matrix_t *m1, const matrix_t *m2);

int mtx_smul2(matrix_t *m_res, const matrix_t *m1, double scalar);
int mtx_sdiv2(matrix_t *m_res, const matrix_t *m1, double scalar);

int mtx_mul(matrix_t *m1_res, const matrix_t *m2_factor);

int mtx_mul2(matrix_t *m_res, const matrix_t *m1, const matrix_t *m2);

int mtx_transpose(matrix_t *m);

int mtx_swap_rows(matrix_t *m, size_t r1, size_t r2);
int mtx_swap_cols(matrix_t *m, size_t c1, size_t c2);

int mtx_scale_row(matrix_t *m, size_t r, double scalar);
int mtx_add_rows(matrix_t *m, size_t r_target, size_t r_source);
int mtx_add_scaled_row(matrix_t *m, size_t r_target, size_t r_source,
                       double scalar);

double mtx_norm(const matrix_t *m);

int mtx_exp(matrix_t *res, const matrix_t *A, double eps);

int mtx_solve_gauss(matrix_t *X_res, const matrix_t *A_orig,
                    const matrix_t *B_orig);

#endif // MATRIX_LIB_H
