#include "matrix_manipulations.h"
#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>

int mtx_transpose(matrix_t *m) {
  if (!m)
    return -1;
  if (m->rows * m->cols == 0)
    return 0;
  if (!m->data)
    return -1;

  matrix_t *temp = mtx_alloc(m->cols, m->rows);
  if (!temp)
    return -1;

  for (size_t i = 0; i < m->rows; ++i) {
    for (size_t j = 0; j < m->cols; ++j) {
      *mtx_ptr(temp, j, i) = *mtx_cptr(m, i, j);
    }
  }

  mtx_move_assign(m, temp);
  return 0;
}

int mtx_swap_rows(matrix_t *m, size_t row1, size_t row2) {
  if (!m || !m->data)
    return -1;
  if (row1 >= m->rows || row2 >= m->rows)
    return -1;

  for (size_t j = 0; j < m->cols; ++j) {
    double temp = *mtx_ptr(m, row1, j);
    *mtx_ptr(m, row1, j) = *mtx_ptr(m, row2, j);
    *mtx_ptr(m, row2, j) = temp;
  }
  return 0;
}

int mtx_swap_cols(matrix_t *m, size_t col1, size_t col2) {
  if (!m || !m->data)
    return -1;
  if (col1 >= m->cols || col2 >= m->cols)
    return -1;

  for (size_t i = 0; i < m->rows; ++i) {
    double temp = *mtx_ptr(m, i, col1);
    *mtx_ptr(m, i, col1) = *mtx_ptr(m, i, col2);
    *mtx_ptr(m, i, col2) = temp;
  }
  return 0;
}

int mtx_scale_row(matrix_t *m, size_t row, double scale) {
  if (!m || !m->data)
    return -1;
  if (row >= m->rows)
    return -1;

  for (size_t j = 0; j < m->cols; ++j) {
    *mtx_ptr(m, row, j) *= scale;
  }
  return 0;
}

int mtx_add_rows(matrix_t *m, size_t dest_row, size_t src_row) {
  if (!m || !m->data)
    return -1;
  if (dest_row >= m->rows || src_row >= m->rows)
    return -1;

  for (size_t j = 0; j < m->cols; ++j) {
    *mtx_ptr(m, dest_row, j) += *mtx_cptr(m, src_row, j);
  }
  return 0;
}

int mtx_add_scaled_row(matrix_t *m, size_t dest_row, size_t src_row, double scale) {
  if (!m || !m->data)
    return -1;
  if (dest_row >= m->rows || src_row >= m->rows)
    return -1;

  for (size_t j = 0; j < m->cols; ++j) {
    *mtx_ptr(m, dest_row, j) += scale * *mtx_cptr(m, src_row, j);
  }
  return 0;
} 