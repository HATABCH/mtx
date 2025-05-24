#ifndef MATRIX_H
#define MATRIX_H

#include <stddef.h>

typedef struct {
  size_t rows;
  size_t cols;
  double *data;
} matrix_t;

matrix_t *mtx_alloc(size_t rows, size_t cols);
matrix_t *mtx_alloc_zero(size_t rows, size_t cols);
matrix_t *mtx_alloc_id(size_t rows, size_t cols);
matrix_t *mtx_copy(const matrix_t *m);
void mtx_free(matrix_t *m);
int mtx_assign(matrix_t *dest, const matrix_t *src);
int mtx_move_assign(matrix_t *dest, matrix_t *src);
void mtx_set_zero(matrix_t *m);
void mtx_set_id(matrix_t *m);
double *mtx_ptr(matrix_t *m, size_t i, size_t j);
const double *mtx_cptr(const matrix_t *m, size_t i, size_t j);
void mtx_print(const matrix_t *m);
void mtx_print_titled(const char *title, const matrix_t *m);

#endif 