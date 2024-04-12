#include "s21_matrix.h"

void support_print_matrix(matrix_t *result) {
  for (int i = 0; i < result->rows; i++) {
    for (int j = 0; j < result->columns; j++) {
      printf("%.lf ", result->matrix[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}

int support_eq_size_matrix(matrix_t *A, matrix_t *B) {
  int res = SUCCESS;
  if ((A->columns != B->columns || A->rows != B->rows)) res = FAILURE;

  return res;
}

void support_search_minor(matrix_t *A, matrix_t *minor, int pass_column,
                          int pass_row) {
  s21_create_matrix(A->rows - 1, A->columns - 1, minor);
  int a = 0, b = 0;

  for (int i = 0; i < A->rows; i++) {
    if (i != pass_row) {
      for (int j = 0; j < A->columns; j++) {
        if (j == pass_column) continue;
        minor->matrix[a][b] = A->matrix[i][j];

        if (b < minor->columns - 1)
          b++;
        else if (a < minor->rows - 1) {
          a++;
          b = 0;
        }
      }
    }
  }
}

int support_be_matrix(matrix_t *A) {
  if (A == NULL || A->matrix == NULL || A->rows <= 0 || A->columns <= 0)
    return FAILURE;
  else
    return SUCCESS;
}