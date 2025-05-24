#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

matrix_t *mtx_alloc(size_t rows, size_t cols) {
  if (rows == 0 || cols == 0)
    return NULL;

  matrix_t *m = (matrix_t *)malloc(sizeof(matrix_t));
  if (!m)
    return NULL;

  m->rows = rows;
  m->cols = cols;
  m->data = (double *)calloc(rows * cols, sizeof(double));
  if (!m->data) {
    free(m);
    return NULL;
  }

  return m;
}

matrix_t *mtx_alloc_zero(size_t rows, size_t cols) {
  return mtx_alloc(rows, cols);
}

matrix_t *mtx_alloc_id(size_t rows, size_t cols) {
  if (rows != cols)
    return NULL;

  matrix_t *m = mtx_alloc(rows, cols);
  if (!m)
    return NULL;

  for (size_t i = 0; i < rows; ++i) {
    *mtx_ptr(m, i, i) = 1.0;
  }

  return m;
}

matrix_t *mtx_copy(const matrix_t *m) {
  if (!m || !m->data)
    return NULL;

  matrix_t *copy = mtx_alloc(m->rows, m->cols);
  if (!copy)
    return NULL;

  memcpy(copy->data, m->data, m->rows * m->cols * sizeof(double));
  return copy;
}

void mtx_free(matrix_t *m) {
  if (m) {
    free(m->data);
    free(m);
  }
}

int mtx_assign(matrix_t *dest, const matrix_t *src) {
  if (!dest || !src || !src->data)
    return -1;
  if (dest->rows != src->rows || dest->cols != src->cols)
    return -1;

  memcpy(dest->data, src->data, src->rows * src->cols * sizeof(double));
  return 0;
}

int mtx_move_assign(matrix_t *dest, matrix_t *src) {
  if (!dest || !src)
    return -1;

  free(dest->data);
  dest->data = src->data;
  dest->rows = src->rows;
  dest->cols = src->cols;

  src->data = NULL;
  src->rows = 0;
  src->cols = 0;

  mtx_free(src);
  return 0;
}

void mtx_set_zero(matrix_t *m) {
  if (m && m->data) {
    memset(m->data, 0, m->rows * m->cols * sizeof(double));
  }
}

void mtx_set_id(matrix_t *m) {
  if (!m || !m->data || m->rows != m->cols)
    return;

  mtx_set_zero(m);
  for (size_t i = 0; i < m->rows; ++i) {
    *mtx_ptr(m, i, i) = 1.0;
  }
}

double *mtx_ptr(matrix_t *m, size_t i, size_t j) {
  if (!m || !m->data || i >= m->rows || j >= m->cols)
    return NULL;
  return &m->data[i * m->cols + j];
}

const double *mtx_cptr(const matrix_t *m, size_t i, size_t j) {
  if (!m || !m->data || i >= m->rows || j >= m->cols)
    return NULL;
  return &m->data[i * m->cols + j];
}

void mtx_print(const matrix_t *m) {
  if (!m || !m->data) {
    printf("NULL matrix\n");
    return;
  }

  for (size_t i = 0; i < m->rows; ++i) {
    for (size_t j = 0; j < m->cols; ++j) {
      printf("%8.3f ", *mtx_cptr(m, i, j));
    }
    printf("\n");
  }
}

void mtx_print_titled(const char *title, const matrix_t *m) {
  printf("%s:\n", title);
  mtx_print(m);
  printf("\n");
} 