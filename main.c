#include "matrix.h"
#include "matrix_manipulations.h"
#include "matrix_operations.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void check_gauss_solution(const matrix_t *A, const matrix_t *X,
                          const matrix_t *B, const char *title);
void test_matrix_exp(void);

void check_gauss_solution(const matrix_t *A, const matrix_t *X,
                          const matrix_t *B, const char *title) {
  printf("--- Проверка решения Гаусса для: %s ---\n", title);

  if (!A || !X || !B) {
    printf("Ошибка: Передан указатель на NULL матрицу для проверки Гаусса.\n");
    return;
  }

  if (A->rows != A->cols || A->rows != X->rows || X->cols != 1 ||
      B->rows != A->rows || B->cols != 1) {
    printf("Ошибка: Несоответствие размерностей матриц для проверки Гаусса.\n");
    return;
  }

  matrix_t *AX = mtx_alloc(A->rows, 1);
  if (!AX) {
    printf("Ошибка: Не удалось выделить память для проверки решения.\n");
    return;
  }

  matrix_t *A_copy = mtx_copy(A);
  if (!A_copy) {
    printf("Ошибка: Не удалось скопировать матрицу A.\n");
    mtx_free(AX);
    return;
  }

  if (mtx_mul(A_copy, X) != 0) {
    printf("Ошибка: Не удалось выполнить умножение A * X.\n");
    mtx_free(AX);
    mtx_free(A_copy);
    return;
  }

  mtx_move_assign(AX, A_copy);

  printf("A * X =\n");
  mtx_print(AX);
  printf("\nB =\n");
  mtx_print(B);

  double diff_norm = 0.0;
  for (size_t i = 0; i < AX->rows; ++i) {
    double diff = *mtx_cptr(AX, i, 0) - *mtx_cptr(B, i, 0);
    diff_norm += diff * diff;
  }
  diff_norm = sqrt(diff_norm);

  printf("Норма разности: %e\n", diff_norm);
  printf("--- Конец проверки ---\n\n");

  mtx_free(AX);
}

void test_matrix_exp(void) {
  printf("\n--- Тест матричной экспоненты ---\n");

  matrix_t *A = mtx_alloc(3, 3);
  if (!A) {
    printf("Ошибка выделения памяти для матрицы A\n");
    return;
  }

  mtx_set_zero(A);
  *mtx_ptr(A, 0, 0) = 1.0;
  *mtx_ptr(A, 1, 1) = 2.0;
  *mtx_ptr(A, 2, 2) = -1.0;

  printf("Исходная матрица A:\n");
  mtx_print(A);

  matrix_t *exp_A = mtx_alloc(3, 3);
  if (!exp_A) {
    printf("Ошибка выделения памяти для exp(A)\n");
    mtx_free(A);
    return;
  }

  if (mtx_exp(A, exp_A) != 0) {
    printf("Ошибка при вычислении матричной экспоненты\n");
    mtx_free(A);
    mtx_free(exp_A);
    return;
  }

  printf("\nРезультат exp(A):\n");
  mtx_print(exp_A);

  mtx_free(A);
  mtx_free(exp_A);
}

int main(void) {
  matrix_t *m1 = mtx_alloc(3, 3);
  matrix_t *m2 = mtx_alloc(3, 3);
  matrix_t *m_copy = NULL;
  matrix_t *m_zero = NULL;
  matrix_t *m_id = NULL;
  matrix_t *gauss_A = NULL;
  matrix_t *gauss_B = NULL;
  matrix_t *gauss_X = NULL;
  matrix_t *aug = NULL;

  if (!m1 || !m2) {
    printf("Ошибка: Не удалось выделить память для матриц.\n");
    goto cleanup;
  }

  double m1_data[] = {2, 1, -1, -3, -1, 2, -2, 1, 2};
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      *mtx_ptr(m1, i, j) = m1_data[i * 3 + j];
    }
  }

  double m2_data[] = {1, 2, 3, 2, 5, 2, 6, -3, 1};
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      *mtx_ptr(m2, i, j) = m2_data[i * 3 + j];
    }
  }

  mtx_print_titled("Матрица a", m1);
  mtx_print_titled("Матрица b", m2);

  m_copy = mtx_copy(m1);
  if (!m_copy) {
    printf("Ошибка: Не удалось скопировать матрицу.\n");
    goto cleanup;
  }
  mtx_print_titled("m_copy (из m1)", m_copy);

  m_id = mtx_alloc_id(3, 3);
  if (!m_id) {
    printf("Ошибка: Не удалось создать единичную матрицу.\n");
    goto cleanup;
  }
  mtx_print_titled("m_id", m_id);

  m_zero = mtx_alloc_zero(3, 3);
  if (!m_zero) {
    printf("Ошибка: Не удалось создать нулевую матрицу.\n");
    goto cleanup;
  }
  mtx_print_titled("m_zero", m_zero);

  if (mtx_assign(m_zero, m_id) != 0) {
    printf("Ошибка: Не удалось присвоить m_zero = m_id.\n");
    goto cleanup;
  }
  mtx_print_titled("m_zero после присваивания из m_id", m_zero);

  if (mtx_add(m1, m2) != 0) {
    printf("Ошибка: Не удалось выполнить сложение матриц.\n");
    goto cleanup;
  }
  mtx_print_titled("a + b", m1);

  if (mtx_sub(m1, m2) != 0) {
    printf("Ошибка: Не удалось выполнить вычитание матриц.\n");
    goto cleanup;
  }
  mtx_print_titled("a - b", m1);

  if (mtx_mul(m1, m2) != 0) {
    printf("Ошибка: Не удалось выполнить умножение матриц.\n");
    goto cleanup;
  }
  mtx_print_titled("a * b", m1);

  if (mtx_div(m1, m2) != 0) {
    printf("Ошибка: Не удалось выполнить деление матриц.\n");
    goto cleanup;
  }
  mtx_print_titled("a / b", m1);

  if (mtx_scale(m1, 2.0) != 0) {
    printf("Ошибка: Не удалось выполнить масштабирование матрицы.\n");
    goto cleanup;
  }
  mtx_print_titled("2 * a", m1);

  if (mtx_transpose(m1) != 0) {
    printf("Ошибка: Не удалось выполнить транспонирование матрицы.\n");
    goto cleanup;
  }
  mtx_print_titled("a^T", m1);

  gauss_A = mtx_alloc(3, 3);
  gauss_B = mtx_alloc(3, 1);
  gauss_X = mtx_alloc(3, 1);

  if (!gauss_A || !gauss_B || !gauss_X) {
    printf("Ошибка: Не удалось выделить память для системы Гаусса.\n");
    goto cleanup;
  }

  double gauss_A1_data[] = {2, 1, -1, -3, -1, 2, -2, 1, 2};
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      *mtx_ptr(gauss_A, i, j) = gauss_A1_data[i * 3 + j];
    }
  }

  double gauss_B1_data[] = {8, -11, -3};
  for (size_t i = 0; i < 3; ++i) {
    *mtx_ptr(gauss_B, i, 0) = gauss_B1_data[i];
  }

  aug = mtx_alloc(3, 4);
  if (!aug) {
    printf("Ошибка: Не удалось выделить память для расширенной матрицы.\n");
    goto cleanup;
  }

  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      *mtx_ptr(aug, i, j) = *mtx_cptr(gauss_A, i, j);
    }
    *mtx_ptr(aug, i, 3) = *mtx_cptr(gauss_B, i, 0);
  }

  if (mtx_gauss_elimination(aug) != 0) {
    printf("Ошибка: Не удалось выполнить метод Гаусса.\n");
    mtx_free(aug);
    aug = NULL;
    goto cleanup;
  }

  for (size_t i = 0; i < 3; ++i) {
    *mtx_ptr(gauss_X, i, 0) = *mtx_cptr(aug, i, 3);
  }

  mtx_print_titled("Гаусс A (Система 1)", gauss_A);
  mtx_print_titled("Гаусс B (Система 1)", gauss_B);
  mtx_print_titled("Гаусс X (Система 1)", gauss_X);

  check_gauss_solution(gauss_A, gauss_X, gauss_B, "Система 1");

  mtx_free(aug);
  aug = NULL;

  double gauss_A2_data[] = {1, 2, 3, 2, 5, 2, 6, -3, 1};
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      *mtx_ptr(gauss_A, i, j) = gauss_A2_data[i * 3 + j];
    }
  }

  double gauss_B2_data[] = {6, 4, 2};
  for (size_t i = 0; i < 3; ++i) {
    *mtx_ptr(gauss_B, i, 0) = gauss_B2_data[i];
  }

  aug = mtx_alloc(3, 4);
  if (!aug) {
    printf("Ошибка: Не удалось выделить память для расширенной матрицы.\n");
    goto cleanup;
  }

  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      *mtx_ptr(aug, i, j) = *mtx_cptr(gauss_A, i, j);
    }
    *mtx_ptr(aug, i, 3) = *mtx_cptr(gauss_B, i, 0);
  }

  if (mtx_gauss_elimination(aug) != 0) {
    printf("Ошибка: Не удалось выполнить метод Гаусса.\n");
    mtx_free(aug);
    aug = NULL;
    goto cleanup;
  }

  for (size_t i = 0; i < 3; ++i) {
    *mtx_ptr(gauss_X, i, 0) = *mtx_cptr(aug, i, 3);
  }

  mtx_print_titled("Гаусс A (Система 2)", gauss_A);
  mtx_print_titled("Гаусс B (Система 2)", gauss_B);
  mtx_print_titled("Гаусс X (Система 2)", gauss_X);

  check_gauss_solution(gauss_A, gauss_X, gauss_B, "Система 2");

  mtx_free(aug);
  aug = NULL;

  printf("\n--- Тест Гаусса ---\n");
  matrix_t *perm_matrix = mtx_alloc(3, 4);
  if (!perm_matrix) {
    printf("Ошибка выделения памяти\n");
    goto cleanup;
  }

  *mtx_ptr(perm_matrix, 0, 0) = 0;
  *mtx_ptr(perm_matrix, 0, 1) = 0;
  *mtx_ptr(perm_matrix, 0, 2) = 1;
  *mtx_ptr(perm_matrix, 0, 3) = 1;
  *mtx_ptr(perm_matrix, 1, 0) = 0;
  *mtx_ptr(perm_matrix, 1, 1) = 1;
  *mtx_ptr(perm_matrix, 1, 2) = 0;
  *mtx_ptr(perm_matrix, 1, 3) = 2;
  *mtx_ptr(perm_matrix, 2, 0) = 1;
  *mtx_ptr(perm_matrix, 2, 1) = 0;
  *mtx_ptr(perm_matrix, 2, 2) = 0;
  *mtx_ptr(perm_matrix, 2, 3) = 3;

  printf("Исходная матрица:\n");
  mtx_print(perm_matrix);

  if (mtx_gauss_elimination(perm_matrix) != 0) {
    printf("Ошибка при выполнении метода Гаусса\n");
    mtx_free(perm_matrix);
    goto cleanup;
  }

  printf("\nМатрица после метода Гаусса:\n");
  mtx_print(perm_matrix);

  printf("\nРешение:\n");
  printf("x1 = %.6f\n", *mtx_cptr(perm_matrix, 0, 3));
  printf("x2 = %.6f\n", *mtx_cptr(perm_matrix, 1, 3));
  printf("x3 = %.6f\n", *mtx_cptr(perm_matrix, 2, 3));

  mtx_free(perm_matrix);

  test_matrix_exp();

cleanup:
  mtx_free(m1);
  mtx_free(m2);
  mtx_free(m_copy);
  mtx_free(m_zero);
  mtx_free(m_id);
  mtx_free(gauss_A);
  mtx_free(gauss_B);
  mtx_free(gauss_X);
  mtx_free(aug);

  return 0;
}
