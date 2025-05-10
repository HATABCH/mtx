#ifndef MATRIX_MANIPULATIONS_H
#define MATRIX_MANIPULATIONS_H

#include "matrix.h"

int mtx_transpose(matrix_t *m);

int mtx_swap_rows(matrix_t *m, size_t row1, size_t row2);
int mtx_swap_cols(matrix_t *m, size_t col1, size_t col2);
int mtx_scale_row(matrix_t *m, size_t row, double scale);
int mtx_add_rows(matrix_t *m, size_t dest_row, size_t src_row);
int mtx_add_scaled_row(matrix_t *m, size_t dest_row, size_t src_row,
                       double scale);

#endif // MATRIX_MANIPULATIONS_H
