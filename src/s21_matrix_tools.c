#include "s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int res = OK;

  if (result == NULL || rows <= 0 || columns <= 0) {
    res = ERROR;
  } else {
    result->rows = rows;
    result->columns = columns;

    result->matrix = calloc(
        result->rows,
        sizeof(double *));  // Даю память на указатели (sizeof(double*)) для
                            // работы с ними в количестве rows (result->rows)

    for (int i = 0; i < result->rows;
         i++)  // даю понять, что в каждой строке будет столько - то столбцов
    {
      result->matrix[i] = calloc(result->columns, sizeof(double));
    }
  }

  return res;
}

void s21_remove_matrix(matrix_t *matrix) {
  for (int i = 0; i < matrix->rows; i++) {
    free(matrix->matrix[i]);
  }
  matrix->rows = 0, matrix->columns = 0;
  free(matrix->matrix);
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int res = SUCCESS;

  if (!(support_eq_size_matrix(A, B)) || (A->columns == 0 && B->rows == 0))
    return FAILURE;
  else
    for (int i = 0; A->rows > i && res; i++) {
      for (int j = 0; A->columns > j && res; j++) {
        if (fabs(A->matrix[i][j] - B->matrix[i][j]) > 1e-07) res = FAILURE;
      }
    }

  return res;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int res = OK;

  if ((!support_be_matrix(A) || !support_be_matrix(B)))
    return ERROR;
  else if (!(support_eq_size_matrix(A, B)))
    return ERROR_CALC;
  else {
    s21_create_matrix(A->rows, A->columns, result);
    for (int i = 0; A->rows > i; i++) {
      for (int j = 0; A->columns > j; j++) {
        result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
      }
    }
  }

  return res;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int res = OK;

  if ((!support_be_matrix(A) || !support_be_matrix(B)))
    return ERROR;
  else if (!(support_eq_size_matrix(A, B)))
    return ERROR_CALC;
  else {
    s21_create_matrix(A->rows, A->columns, result);
    for (int i = 0; A->rows > i; i++) {
      for (int j = 0; A->columns > j; j++) {
        result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
      }
    }
  }

  return res;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int res = OK;

  if (!support_be_matrix(A))
    return ERROR;
  else {
    s21_create_matrix(A->columns, A->rows, result);
    for (int i = 0; A->rows > i; i++) {
      for (int j = 0; A->columns > j; j++) {
        result->matrix[j][i] = A->matrix[i][j];
      }
    }
  }

  return res;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int res = OK;

  if (!support_be_matrix(A))
    return ERROR;
  else {
    s21_create_matrix(A->rows, A->columns, result);
    for (int i = 0; A->rows > i; i++) {
      for (int j = 0; A->columns > j; j++) {
        result->matrix[i][j] = A->matrix[i][j] * number;
      }
    }
  }

  return res;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int res = OK;

  if (!support_be_matrix(A) || !support_be_matrix(B)) {
    return ERROR;
  }
  if (A->columns == B->rows) {
    s21_create_matrix(A->rows, B->columns, result);
    for (int i = 0; i < result->rows; i++) {
      for (int j = 0; j < result->columns; j++) {
        for (int k = 0; k < A->columns; k++) {
          result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
        }
      }
    }
  } else {
    res = ERROR_CALC;
  }
  return res;
}

int s21_determinant(matrix_t *A, double *result) {
  *result = 0.0;
  int res = OK;
  // matrix_t *minor = calloc(1, sizeof(matrix_t));

  if (A->columns != A->rows)
    return ERROR_CALC;
  else if (!support_be_matrix(A))
    return ERROR;
  else if (A->rows == 1 && A->columns == 1)
    *result = A->matrix[0][0];
  else if (A->rows == 2 && A->columns == 2)
    *result =
        A->matrix[0][0] * A->matrix[1][1] - A->matrix[1][0] * A->matrix[0][1];
  if (A->rows > 2) {
    int sign = 1;
    for (int i = 0; A->rows > i; i++) {
      matrix_t minor;
      // matrix_t *minor = calloc(1, sizeof(matrix_t));
      double result_two = 0.0;
      support_search_minor(A, &minor, i, 0);
      s21_determinant(&minor, &result_two);
      *result += sign * A->matrix[0][i] * result_two;
      sign *= -1;
      s21_remove_matrix(&minor);
    }
  }
  // s21_remove_matrix(minor);
  // free(minor);
  return res;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int res = OK;

  if (support_be_matrix(A) == 0)
    return ERROR;
  else if (A->rows != A->columns)
    return ERROR_CALC;

  if ((s21_create_matrix(A->rows, A->columns, result) == 0)) {
    // matrix_t *minor = calloc(1, sizeof(matrix_t));
    // s21_create_matrix(A->rows, A->columns, minor);

    for (int i = 0; A->rows > i; i++) {
      for (int j = 0; A->columns > j; j++) {
        double result_det = 0.0;
        matrix_t *minor = calloc(1, sizeof(matrix_t));
        support_search_minor(A, minor, i, j);
        s21_determinant(minor, &result_det);
        result->matrix[j][i] = pow(-1, i + j) * result_det;
        s21_remove_matrix(minor);
        free(minor);
      }
    }
    // free(minor);
    // s21_remove_matrix(minor);
  }
  // s21_remove_matrix(result);
  // free(result);
  return res;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int res = OK;
  double result_det = 0.0;

  if (!support_be_matrix(A)) return ERROR;
  if (A->columns != A->rows)
    return ERROR_CALC;
  else {
    s21_determinant(A, &result_det);
    if (result_det == 0) res = ERROR_CALC;

    matrix_t calc;       // = calloc(1, sizeof(matrix_t));
    matrix_t transpose;  //= calloc(1, sizeof(matrix_t));
                         // s21_create_matrix(A->rows, A->columns, &calc);
    // s21_create_matrix(A->rows, A->columns, result);
    // s21_create_matrix(A->rows, A->columns, &transpose);

    // for (int i = 0; A->rows > i; i++) {
    // for (int j = 0; A->columns > j; j++) {
    s21_calc_complements(A, &calc);
    s21_transpose(&calc, &transpose);
    s21_mult_number(&transpose, 1 / result_det, result);
    //}
    //}

    s21_remove_matrix(&calc);
    s21_remove_matrix(&transpose);
    // free(calc);
    // free(transpose);
  }
  return res;
}
