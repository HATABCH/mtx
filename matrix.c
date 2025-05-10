#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

matrix_t *mtx_alloc(size_t h, size_t w) {
  matrix_t *m = (matrix_t *)malloc(sizeof(matrix_t));
  if (!m) {
    return NULL;
  }

  m->rows = h;
  m->cols = w;
  if (h * w == 0) {
    m->data = NULL;
  } else {
    m->data = (double *)malloc(h * w * sizeof(double));
    if (!m->data) {
      free(m);
      return NULL;
    }
  }
  return m;
}

matrix_t *mtx_copy(const matrix_t *src) {
  if (!src) {
    return NULL;
  }
  matrix_t *dst = mtx_alloc(src->rows, src->cols);
  if (!dst) {
    return NULL;
  }
  if (src->data && dst->data) {
    memcpy(dst->data, src->data, src->rows * src->cols * sizeof(double));
  }
  return dst;
}

void mtx_free(matrix_t *m) {
  if (m) {
    free(m->data);
    free(m);
  }
}

double *mtx_ptr(matrix_t *m, size_t r, size_t c) {
  return &m->data[r * m->cols + c];
}

const double *mtx_cptr(const matrix_t *m, size_t r, size_t c) {
  return &m->data[r * m->cols + c];
}

size_t mtx_get_rows(const matrix_t *m) { return m ? m->rows : 0; }

size_t mtx_get_cols(const matrix_t *m) { return m ? m->cols : 0; }

void mtx_set_zero(matrix_t *m) {
  if (m && m->data) {
    memset(m->data, 0, m->rows * m->cols * sizeof(double));
  }
}

void mtx_set_id(matrix_t *m) {
  if (!m)
    return;
  mtx_set_zero(m);
  if (m->data) {
    size_t min_dim = (m->rows < m->cols) ? m->rows : m->cols;
    for (size_t i = 0; i < min_dim; ++i) {
      *mtx_ptr(m, i, i) = 1.0;
    }
  }
}

matrix_t *mtx_alloc_zero(size_t h, size_t w) {
  matrix_t *m = mtx_alloc(h, w);
  if (m) {
    mtx_set_zero(m);
  }
  return m;
}

matrix_t *mtx_alloc_id(size_t h, size_t w) {
  matrix_t *m = mtx_alloc(h, w);
  if (m) {
    mtx_set_id(m);
  }
  return m;
}

int mtx_assign(matrix_t *m1_dst, const matrix_t *m2_src) {
  if (!m1_dst || !m2_src)
    return -1;
  if (m1_dst->rows != m2_src->rows || m1_dst->cols != m2_src->cols) {
    return -1;
  }
  if (m1_dst->rows * m1_dst->cols == 0) {
    return 0;
  }
  if (!m1_dst->data || !m2_src->data) {
    fprintf(stderr, "Предупреждение: mtx_assign вызван для матрицы ненулевого размера с NULL данными.\n");
    return -1;
  }
  memcpy(m1_dst->data, m2_src->data,
         m1_dst->rows * m1_dst->cols * sizeof(double));
  return 0;
}

void mtx_move_assign(matrix_t *target, matrix_t *source_to_consume) {
  if (!target || !source_to_consume)
    return;
  free(target->data);
  target->data = source_to_consume->data;
  target->rows = source_to_consume->rows;
  target->cols = source_to_consume->cols;
  source_to_consume->data = NULL;
  source_to_consume->rows = 0;
  source_to_consume->cols = 0;
}

int mtx_read(matrix_t *m) {
  if (!m)
    return -1;
  printf("Введите матрицу размером %zu x %zu (построчно):\n", mtx_get_rows(m),
         mtx_get_cols(m));
  if (mtx_get_rows(m) * mtx_get_cols(m) > 0 && !m->data) {
    fprintf(stderr, "Ошибка: mtx_read вызван для непустой матрицы с NULL буфером данных.\n");
    return -1;
  }
  for (size_t r = 0; r < mtx_get_rows(m); ++r) {
    for (size_t c = 0; c < mtx_get_cols(m); ++c) {
      if (scanf("%lf", mtx_ptr(m, r, c)) != 1) {
        fprintf(stderr, "Ошибка чтения элемента [%zu][%zu]\n", r, c);
        int ch;
        while ((ch = getchar()) != '\n' && ch != EOF)
          ;
        return -1;
      }
    }
  }
  return 0;
}

void mtx_print(const matrix_t *m) {
  if (!m) {
    printf("(указатель на матрицу равен NULL)\n");
    return;
  }
  if (mtx_get_rows(m) * mtx_get_cols(m) == 0) {
    printf("Матрица %zu x %zu (пустая или данные равны NULL)\n", mtx_get_rows(m),
           mtx_get_cols(m));
    return;
  }
  if (!m->data) {
    printf("Матрица %zu x %zu (данные неожиданно равны NULL для непустой матрицы)\n",
        mtx_get_rows(m), mtx_get_cols(m));
    return;
  }
  for (size_t r = 0; r < mtx_get_rows(m); ++r) {
    for (size_t c = 0; c < mtx_get_cols(m); ++c) {
      printf("%8.3f ", *mtx_cptr(m, r, c));
    }
    printf("\n");
  }
}

void mtx_print_titled(const char *title, const matrix_t *m) {
  printf("%s:\n", title);
  mtx_print(m);
  printf("\n");
} 