#include "matrix_operations.h"
#include "matrix.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int mtx_add(matrix_t *m1, const matrix_t *m2) {
  if (!m1 || !m2)
    return -1;
  if (m1->rows != m2->rows || m1->cols != m2->cols)
    return -1;
  if (m1->rows * m1->cols == 0)
    return 0;
  if (!m1->data || !m2->data)
    return -1;

  for (size_t i = 0; i < m1->rows * m1->cols; ++i) {
    m1->data[i] += m2->data[i];
  }
  return 0;
}

int mtx_sub(matrix_t *m1, const matrix_t *m2) {
  if (!m1 || !m2)
    return -1;
  if (m1->rows != m2->rows || m1->cols != m2->cols)
    return -1;
  if (m1->rows * m1->cols == 0)
    return 0;
  if (!m1->data || !m2->data)
    return -1;

  for (size_t i = 0; i < m1->rows * m1->cols; ++i) {
    m1->data[i] -= m2->data[i];
  }
  return 0;
}

int mtx_mul(matrix_t *m1, const matrix_t *m2) {
  if (!m1 || !m2)
    return -1;
  if (m1->cols != m2->rows)
    return -1;
  if (m1->rows * m1->cols == 0 || m2->rows * m2->cols == 0)
    return 0;
  if (!m1->data || !m2->data)
    return -1;

  matrix_t *temp = mtx_alloc(m1->rows, m2->cols);
  if (!temp)
    return -1;

  for (size_t i = 0; i < m1->rows; ++i) {
    for (size_t j = 0; j < m2->cols; ++j) {
      double sum = 0.0;
      for (size_t k = 0; k < m1->cols; ++k) {
        sum += *mtx_cptr(m1, i, k) * *mtx_cptr(m2, k, j);
      }
      *mtx_ptr(temp, i, j) = sum;
    }
  }

  mtx_move_assign(m1, temp);
  return 0;
}

int mtx_div(matrix_t *m1, const matrix_t *m2) {
  if (!m1 || !m2)
    return -1;
  if (m1->cols != m2->rows || m2->rows != m2->cols)
    return -1;
  if (m1->rows * m1->cols == 0 || m2->rows * m2->cols == 0)
    return 0;
  if (!m1->data || !m2->data)
    return -1;

  matrix_t *inv = mtx_alloc(m2->rows, m2->cols);
  if (!inv)
    return -1;

  if (mtx_inverse(m2, inv) != 0) {
    mtx_free(inv);
    return -1;
  }

  matrix_t *temp = mtx_alloc(m1->rows, m2->cols);
  if (!temp) {
    mtx_free(inv);
    return -1;
  }

  matrix_t *m1_copy = mtx_copy(m1);
  if (!m1_copy) {
    mtx_free(inv);
    mtx_free(temp);
    return -1;
  }

  if (mtx_mul(m1_copy, inv) != 0) {
    mtx_free(inv);
    mtx_free(temp);
    mtx_free(m1_copy);
    return -1;
  }

  mtx_move_assign(temp, m1_copy);
  mtx_move_assign(m1, temp);
  mtx_free(inv);
  return 0;
}

int mtx_add_scaled(matrix_t *m1, const matrix_t *m2, double scale) {
  if (!m1 || !m2)
    return -1;
  if (m1->rows != m2->rows || m1->cols != m2->cols)
    return -1;
  if (m1->rows * m1->cols == 0)
    return 0;
  if (!m1->data || !m2->data)
    return -1;

  for (size_t i = 0; i < m1->rows * m1->cols; ++i) {
    m1->data[i] += scale * m2->data[i];
  }
  return 0;
}

int mtx_scale(matrix_t *m, double scale) {
  if (!m)
    return -1;
  if (m->rows * m->cols == 0)
    return 0;
  if (!m->data)
    return -1;

  for (size_t i = 0; i < m->rows * m->cols; ++i) {
    m->data[i] *= scale;
  }
  return 0;
}

double mtx_norm(const matrix_t *m) {
  if (!m || !m->data)
    return 0.0;

  double sum = 0.0;
  for (size_t i = 0; i < m->rows * m->cols; ++i) {
    sum += m->data[i] * m->data[i];
  }
  return sqrt(sum);
}

int mtx_inverse(const matrix_t *m, matrix_t *inv) {
  if (!m || !inv)
    return -1;
  if (m->rows != m->cols || inv->rows != inv->cols || m->rows != inv->rows)
    return -1;
  if (m->rows * m->cols == 0)
    return 0;
  if (!m->data || !inv->data)
    return -1;

  matrix_t *aug = mtx_alloc(m->rows, 2 * m->cols);
  if (!aug)
    return -1;

  for (size_t i = 0; i < m->rows; ++i) {
    for (size_t j = 0; j < m->cols; ++j) {
      *mtx_ptr(aug, i, j) = *mtx_cptr(m, i, j);
    }
    *mtx_ptr(aug, i, i + m->cols) = 1.0;
  }

  if (mtx_gauss_elimination(aug) != 0) {
    mtx_free(aug);
    return -1;
  }

  for (size_t i = 0; i < m->rows; ++i) {
    for (size_t j = 0; j < m->cols; ++j) {
      *mtx_ptr(inv, i, j) = *mtx_cptr(aug, i, j + m->cols);
    }
  }

  mtx_free(aug);
  return 0;
}

int mtx_gauss_elimination(matrix_t *m) {
  if (!m || !m->data)
    return -1;
  if (m->rows == 0 || m->cols == 0)
    return 0;

  size_t n = m->rows;
  size_t m_cols = m->cols;

  for (size_t i = 0; i < n; ++i) {
    size_t max_row = i;
    double max_val = fabs(*mtx_cptr(m, i, i));

    for (size_t j = i + 1; j < n; ++j) {
      double val = fabs(*mtx_cptr(m, j, i));
      if (val > max_val) {
        max_val = val;
        max_row = j;
      }
    }

    if (max_val < EPSILON) {
      return -1;
    }

    if (max_row != i) {
      for (size_t j = 0; j < m_cols; ++j) {
        double temp = *mtx_ptr(m, i, j);
        *mtx_ptr(m, i, j) = *mtx_ptr(m, max_row, j);
        *mtx_ptr(m, max_row, j) = temp;
      }
    }

    for (size_t j = i + 1; j < n; ++j) {
      double factor = *mtx_cptr(m, j, i) / *mtx_cptr(m, i, i);
      for (size_t k = i; k < m_cols; ++k) {
        *mtx_ptr(m, j, k) -= factor * *mtx_cptr(m, i, k);
      }
    }
  }

  for (int i = n - 1; i >= 0; --i) {
    for (int j = i - 1; j >= 0; --j) {
      double factor = *mtx_cptr(m, j, i) / *mtx_cptr(m, i, i);
      for (size_t k = i; k < m_cols; ++k) {
        *mtx_ptr(m, j, k) -= factor * *mtx_cptr(m, i, k);
      }
    }
  }

  for (size_t i = 0; i < n; ++i) {
    double diag = *mtx_cptr(m, i, i);
    for (size_t j = 0; j < m_cols; ++j) {
      *mtx_ptr(m, i, j) /= diag;
    }
  }

  return 0;
}
