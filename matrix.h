#ifndef MATRIX_H
#define MATRIX_H

#include <stddef.h>

typedef struct {
  double *data;
  size_t rows;
  size_t cols;
} matrix_t;

matrix_t *mtx_alloc(size_t h, size_t w);
matrix_t *mtx_copy(const matrix_t *src);
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

#endif 