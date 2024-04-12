#include <check.h>

#include "s21_matrix.h"

void s21_init_matrix(double number, matrix_t *A) {
  for (int x = 0; x < A->rows; x += 1) {
    for (int y = 0; y < A->columns; number += 1.0, y += 1)
      A->matrix[x][y] = number;
  }
}

START_TEST(s21_create_matrix_1) {
  int res = 0;
  matrix_t A = {0};

  res = s21_create_matrix(1, 1, &A);
  ck_assert_int_eq(res, OK);
  s21_remove_matrix(&A);
}

START_TEST(s21_create_matrix_2) {
  int res = 0;
  matrix_t A = {0};

  res = s21_create_matrix(5, 5, &A);
  ck_assert_int_eq(res, OK);
  s21_remove_matrix(&A);
}

START_TEST(s21_create_matrix_3) {
  int res = 0;
  matrix_t A = {0};

  res = s21_create_matrix(5, -1, &A);
  ck_assert_int_eq(res, ERROR);
  s21_remove_matrix(&A);
}

START_TEST(s21_create_matrix_4) {
  int res = 0;
  matrix_t A = {0};

  res = s21_create_matrix(-5, 5, &A);
  ck_assert_int_eq(res, ERROR);
  s21_remove_matrix(&A);
}

START_TEST(s21_create_matrix_5) {
  int res = 0;
  matrix_t A = {0};

  res = s21_create_matrix(0, 5, &A);
  ck_assert_int_eq(res, ERROR);
  s21_remove_matrix(&A);
}

START_TEST(create_1) {
  const int rows = rand() % 100 + 1;
  const int cols = rand() % 100 + 1;
  matrix_t m = {0};
  s21_create_matrix(rows, cols, &m);
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      m.matrix[i][j] = 0;
      ck_assert_ldouble_eq_tol(0, m.matrix[i][j], 1e-07);
    }
  }
  ck_assert_int_eq(m.rows, rows);
  ck_assert_int_eq(m.columns, cols);
  s21_remove_matrix(&m);
}
END_TEST

START_TEST(create_2) {
  const int rows = 0;
  const int cols = 10;
  matrix_t m = {0};
  ck_assert_int_eq(s21_create_matrix(rows, cols, &m), ERROR);
}
END_TEST

START_TEST(create_3) {
  const int rows = 10;
  const int cols = 0;

  matrix_t m = {0};
  ck_assert_int_eq(s21_create_matrix(rows, cols, &m), ERROR);
}
END_TEST

START_TEST(s21_eq_matrix_1) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};

  s21_create_matrix(3, 2, &A);
  s21_init_matrix(1.0, &A);

  s21_create_matrix(3, 2, &B);
  s21_init_matrix(1.0, &B);

  res = s21_eq_matrix(&A, &B);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_eq_matrix_2) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};

  s21_create_matrix(3, 2, &A);
  s21_init_matrix(1.0, &A);

  s21_create_matrix(3, 2, &B);
  s21_init_matrix(2.0, &B);

  res = s21_eq_matrix(&A, &B);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);

  ck_assert_int_eq(res, FAILURE);
}
END_TEST

START_TEST(s21_eq_matrix_3) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};

  s21_create_matrix(4, 2, &A);
  s21_init_matrix(1.0, &A);

  s21_create_matrix(3, 2, &B);
  s21_init_matrix(1.0, &B);

  res = s21_eq_matrix(&A, &B);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);

  ck_assert_int_eq(res, FAILURE);
}
END_TEST

START_TEST(s21_eq_matrix_4) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};

  s21_create_matrix(3, 2, &A);
  s21_init_matrix(1.0, &A);

  s21_create_matrix(3, 2, &B);
  s21_init_matrix(1.0, &B);

  A.matrix[0][0] = 2.0;
  A.matrix[0][1] = 4.0;
  A.matrix[1][0] = 6.0;
  A.matrix[1][1] = 12.0;
  A.matrix[2][0] = 18.0;
  A.matrix[2][1] = 20.0;

  B.matrix[0][0] = 2.0;
  B.matrix[0][1] = 4.0;
  B.matrix[1][0] = 6.0;
  B.matrix[1][1] = 12.0;
  B.matrix[2][0] = 18.0;
  B.matrix[2][1] = 20.000000001;

  res = s21_eq_matrix(&A, &B);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_eq_matrix_5) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};

  s21_create_matrix(3, 2, &A);
  s21_init_matrix(1.0, &A);

  s21_create_matrix(3, 2, &B);
  s21_init_matrix(1.0, &B);

  A.matrix[0][0] = 2.0;
  A.matrix[0][1] = 4.0;
  A.matrix[1][0] = 6.0;
  A.matrix[1][1] = 12.0;
  A.matrix[2][0] = 18.0;
  A.matrix[2][1] = 20.0;

  B.matrix[0][0] = 2.0;
  B.matrix[0][1] = 4.0;
  B.matrix[1][0] = 6.0;
  B.matrix[1][1] = 12.0;
  B.matrix[2][0] = 18.0;
  B.matrix[2][1] = 20.0000001;

  res = s21_eq_matrix(&A, &B);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);

  ck_assert_int_eq(res, FAILURE);
}
END_TEST

START_TEST(zero_matrix) {
  matrix_t A = {0};
  matrix_t B = {0};
  int result = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(0, result);
}

START_TEST(zero_matrix_1) {
  matrix_t A = {0};
  matrix_t B = {0};
  s21_create_matrix(0, 0, &A);
  s21_create_matrix(0, 0, &B);
  int result = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(0, result);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}

START_TEST(casual_matrix_1) {
  matrix_t A = {0};
  matrix_t B = {0};
  s21_create_matrix(1, 1, &A);
  s21_create_matrix(1, 1, &B);
  A.matrix[0][0] = 1;
  B.matrix[0][0] = 1;
  int result = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(1, result);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}

START_TEST(casual_matrix_2) {
  matrix_t A = {0};
  matrix_t B = {0};
  s21_create_matrix(1, 1, &A);
  s21_create_matrix(1, 1, &B);
  A.matrix[0][0] = 1;
  B.matrix[0][0] = 2;
  int result = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(0, result);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}

START_TEST(casual_matrix_3) {
  matrix_t A = {0};
  matrix_t B = {0};
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);
  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[1][0] = 3;
  A.matrix[1][1] = 4;
  B.matrix[0][0] = 1;
  B.matrix[0][1] = 2;
  B.matrix[1][0] = 3;
  B.matrix[1][1] = 4;
  int result = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(1, result);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}

START_TEST(casual_matrix_4) {
  matrix_t A = {0};
  matrix_t B = {0};
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);
  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[1][0] = 3;
  A.matrix[1][1] = 4;
  B.matrix[0][0] = 1;
  B.matrix[0][1] = 2;
  B.matrix[1][0] = 3;
  B.matrix[1][1] = 4;
  int result = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(1, result);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}

START_TEST(casual_matrix_5) {
  matrix_t A = {0};
  matrix_t B = {0};
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);
  A.matrix[0][0] = 1.01;
  A.matrix[0][1] = 2;
  A.matrix[1][0] = 3.05;
  A.matrix[1][1] = 4;
  B.matrix[0][0] = 1.01;
  B.matrix[0][1] = 2;
  B.matrix[1][0] = 3.05;
  B.matrix[1][1] = 4;
  int result = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(1, result);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}

START_TEST(casual_matrix_6) {
  matrix_t A = {0};
  matrix_t B = {0};
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);
  A.matrix[0][0] = 1.01;
  A.matrix[0][1] = -2;
  A.matrix[1][0] = 3.05;
  A.matrix[1][1] = -4;
  B.matrix[0][0] = 1.01;
  B.matrix[0][1] = -2;
  B.matrix[1][0] = 3.05;
  B.matrix[1][1] = -4;
  int result = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(1, result);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}

START_TEST(casual_matrix_7) {
  matrix_t A = {0};
  matrix_t B = {0};
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);
  A.matrix[0][0] = 1.00000000234;
  A.matrix[0][1] = 2;
  A.matrix[1][0] = 3.05;
  A.matrix[1][1] = 4;
  B.matrix[0][0] = 1;
  B.matrix[0][1] = 2;
  B.matrix[1][0] = 3.05;
  B.matrix[1][1] = 4;
  int result = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(1, result);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}

START_TEST(casual_matrix_8) {
  matrix_t A = {0};
  matrix_t B = {0};
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);
  A.matrix[0][0] = 1.0001;
  A.matrix[0][1] = 2;
  A.matrix[1][0] = 3.05;
  A.matrix[1][1] = 4;
  B.matrix[0][0] = 1;
  B.matrix[0][1] = 2;
  B.matrix[1][0] = 3.05;
  B.matrix[1][1] = 4;
  int result = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(0, result);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}

START_TEST(s21_sum_matrix_1) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(4, 4, &A);
  s21_init_matrix(1.0, &A);

  s21_create_matrix(4, 4, &B);
  s21_init_matrix(1.0, &B);

  res = s21_sum_matrix(&A, &B, &Z);
  ck_assert_int_eq(res, OK);

  s21_create_matrix(4, 4, &X);

  X.matrix[0][0] = 2.0;
  X.matrix[0][1] = 4.0;
  X.matrix[0][2] = 6.0;
  X.matrix[0][3] = 8.0;
  X.matrix[1][0] = 10.0;
  X.matrix[1][1] = 12.0;
  X.matrix[1][2] = 14.0;
  X.matrix[1][3] = 16.0;
  X.matrix[2][0] = 18.0;
  X.matrix[2][1] = 20.0;
  X.matrix[2][2] = 22.0;
  X.matrix[2][3] = 24.0;
  X.matrix[3][0] = 26.0;
  X.matrix[3][1] = 28.0;
  X.matrix[3][2] = 30.0;
  X.matrix[3][3] = 32.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_sum_matrix_2) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t Z = {0};

  s21_create_matrix(5, 4, &A);
  s21_init_matrix(0.0, &A);

  s21_create_matrix(4, 4, &B);
  s21_init_matrix(2.0, &B);

  res = s21_sum_matrix(&A, &B, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, ERROR_CALC);
}
END_TEST

START_TEST(s21_sum_matrix_3) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(4, 4, &A);
  s21_init_matrix(-10.0, &A);

  s21_create_matrix(4, 4, &B);
  s21_init_matrix(5.0, &B);

  res = s21_sum_matrix(&A, &B, &Z);
  ck_assert_int_eq(res, OK);

  s21_create_matrix(4, 4, &X);

  X.matrix[0][0] = -5.0;
  X.matrix[0][1] = -3.0;
  X.matrix[0][2] = -1.0;
  X.matrix[0][3] = 1.0;
  X.matrix[1][0] = 3.0;
  X.matrix[1][1] = 5.0;
  X.matrix[1][2] = 7.0;
  X.matrix[1][3] = 9.0;
  X.matrix[2][0] = 11.0;
  X.matrix[2][1] = 13.0;
  X.matrix[2][2] = 15.0;
  X.matrix[2][3] = 17.0;
  X.matrix[3][0] = 19.0;
  X.matrix[3][1] = 21.0;
  X.matrix[3][2] = 23.0;
  X.matrix[3][3] = 25.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_sum_matrix_4) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t Z = {0};

  s21_create_matrix(2, 4, &A);
  s21_init_matrix(-3.0, &A);

  s21_create_matrix(0, 0, &B);

  res = s21_sum_matrix(&A, &B, &Z);
  ck_assert_int_eq(res, ERROR);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&Z);
}
END_TEST

START_TEST(sum_matrix2) {
  int rows = rand() % 100 + 1;
  int cols = rand() % 100 + 1;
  rows = -rows;
  cols = -cols;
  matrix_t m = {0};
  s21_create_matrix(rows, cols, &m);
  matrix_t mtx = {0};
  s21_create_matrix(rows, cols, &mtx);
  matrix_t check = {0};
  s21_create_matrix(rows, cols, &check);

  matrix_t res = {0};

  ck_assert_int_eq(s21_sum_matrix(&m, &mtx, &res), ERROR);
  s21_remove_matrix(&m);
  s21_remove_matrix(&mtx);
  s21_remove_matrix(&res);
  s21_remove_matrix(&check);
}
END_TEST

START_TEST(s21_sub_matrix_1) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(4, 4, &A);
  s21_init_matrix(1.0, &A);

  s21_create_matrix(4, 4, &B);
  s21_init_matrix(1.0, &B);

  res = s21_sub_matrix(&A, &B, &Z);
  ck_assert_int_eq(res, OK);

  s21_create_matrix(4, 4, &X);

  X.matrix[0][0] = 0.0;
  X.matrix[0][1] = 0.0;
  X.matrix[0][2] = 0.0;
  X.matrix[0][3] = 0.0;
  X.matrix[1][0] = 0.0;
  X.matrix[1][1] = 0.0;
  X.matrix[1][2] = 0.0;
  X.matrix[1][3] = 0.0;
  X.matrix[2][0] = 0.0;
  X.matrix[2][1] = 0.0;
  X.matrix[2][2] = 0.0;
  X.matrix[2][3] = 0.0;
  X.matrix[3][0] = 0.0;
  X.matrix[3][1] = 0.0;
  X.matrix[3][2] = 0.0;
  X.matrix[3][3] = 0.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_sub_matrix_2) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(4, 4, &A);
  s21_init_matrix(0.0, &A);

  s21_create_matrix(4, 4, &B);
  s21_init_matrix(-3.0, &B);

  res = s21_sub_matrix(&A, &B, &Z);
  ck_assert_int_eq(res, OK);

  s21_create_matrix(4, 4, &X);

  X.matrix[0][0] = 3.0;
  X.matrix[0][1] = 3.0;
  X.matrix[0][2] = 3.0;
  X.matrix[0][3] = 3.0;
  X.matrix[1][0] = 3.0;
  X.matrix[1][1] = 3.0;
  X.matrix[1][2] = 3.0;
  X.matrix[1][3] = 3.0;
  X.matrix[2][0] = 3.0;
  X.matrix[2][1] = 3.0;
  X.matrix[2][2] = 3.0;
  X.matrix[2][3] = 3.0;
  X.matrix[3][0] = 3.0;
  X.matrix[3][1] = 3.0;
  X.matrix[3][2] = 3.0;
  X.matrix[3][3] = 3.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_sub_matrix_3) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t Z = {0};

  s21_create_matrix(2, 4, &A);
  s21_init_matrix(-3.0, &A);

  s21_create_matrix(4, 4, &B);
  s21_init_matrix(5.0, &B);

  res = s21_sub_matrix(&A, &B, &Z);
  ck_assert_int_eq(res, ERROR_CALC);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&Z);
}
END_TEST

START_TEST(s21_sub_matrix_4) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t Z = {0};

  s21_create_matrix(2, 4, &A);
  s21_init_matrix(-3.0, &A);

  res = s21_sub_matrix(&A, &B, &Z);
  ck_assert_int_eq(res, ERROR);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&Z);
}
END_TEST

START_TEST(sub_matrix2) {
  int rows = rand() % 100 + 1;
  int cols = rand() % 100 + 1;
  rows = -rows;
  cols = -cols;
  matrix_t m = {0};
  s21_create_matrix(rows, cols, &m);
  matrix_t mtx = {0};
  s21_create_matrix(rows, cols, &mtx);
  matrix_t check = {0};
  s21_create_matrix(rows, cols, &check);

  matrix_t res = {0};
  ck_assert_int_eq(s21_sub_matrix(&m, &mtx, &res), ERROR);

  s21_remove_matrix(&m);
  s21_remove_matrix(&mtx);
  s21_remove_matrix(&res);
  s21_remove_matrix(&check);
}
END_TEST

START_TEST(s21_mult_number_1) {
  int res = 0;
  double number = 3.0;
  matrix_t A = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(4, 4, &A);

  s21_init_matrix(1.0, &A);

  res = s21_mult_number(&A, number, &Z);
  ck_assert_int_eq(res, OK);

  s21_create_matrix(4, 4, &X);

  X.matrix[0][0] = 3.0;
  X.matrix[0][1] = 6.0;
  X.matrix[0][2] = 9.0;
  X.matrix[0][3] = 12.0;
  X.matrix[1][0] = 15.0;
  X.matrix[1][1] = 18.0;
  X.matrix[1][2] = 21.0;
  X.matrix[1][3] = 24.0;
  X.matrix[2][0] = 27.0;
  X.matrix[2][1] = 30.0;
  X.matrix[2][2] = 33.0;
  X.matrix[2][3] = 36.0;
  X.matrix[3][0] = 39.0;
  X.matrix[3][1] = 42.0;
  X.matrix[3][2] = 45.0;
  X.matrix[3][3] = 48.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_mult_number_2) {
  int res = 0;
  double number = 5.0;
  matrix_t A = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(4, 4, &A);
  s21_init_matrix(1.0, &A);

  res = s21_mult_number(&A, number, &Z);
  ck_assert_int_eq(res, OK);

  s21_create_matrix(4, 4, &X);

  X.matrix[0][0] = 5.0;
  X.matrix[0][1] = 10.0;
  X.matrix[0][2] = 15.0;
  X.matrix[0][3] = 20.0;
  X.matrix[1][0] = 25.0;
  X.matrix[1][1] = 30.0;
  X.matrix[1][2] = 35.0;
  X.matrix[1][3] = 40.0;
  X.matrix[2][0] = 45.0;
  X.matrix[2][1] = 50.0;
  X.matrix[2][2] = 55.0;
  X.matrix[2][3] = 60.0;
  X.matrix[3][0] = 65.0;
  X.matrix[3][1] = 70.0;
  X.matrix[3][2] = 75.0;
  X.matrix[3][3] = 80.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_mult_number_3) {
  int res = 0;
  double number = 2.0;
  matrix_t A = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(4, 2, &A);
  s21_init_matrix(2.0, &A);

  res = s21_mult_number(&A, number, &Z);
  ck_assert_int_eq(res, OK);

  s21_create_matrix(4, 2, &X);

  X.matrix[0][0] = 4.0;
  X.matrix[0][1] = 6.0;
  X.matrix[1][0] = 8.0;
  X.matrix[1][1] = 10.0;
  X.matrix[2][0] = 12.0;
  X.matrix[2][1] = 14.0;
  X.matrix[3][0] = 16.0;
  X.matrix[3][1] = 18.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_mult_number_4) {
  int res = 0;
  double number = -2.0;
  matrix_t A = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(4, 2, &A);
  s21_init_matrix(2.0, &A);

  res = s21_mult_number(&A, number, &Z);
  ck_assert_int_eq(res, OK);

  s21_create_matrix(4, 2, &X);

  X.matrix[0][0] = -4.0;
  X.matrix[0][1] = -6.0;
  X.matrix[1][0] = -8.0;
  X.matrix[1][1] = -10.0;
  X.matrix[2][0] = -12.0;
  X.matrix[2][1] = -14.0;
  X.matrix[3][0] = -16.0;
  X.matrix[3][1] = -18.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_mult_number_5) {
  int res = 0;
  double number = 2.0;
  matrix_t A = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(1, 1, &A);
  s21_init_matrix(2.0, &A);

  res = s21_mult_number(&A, number, &Z);
  ck_assert_int_eq(res, OK);

  s21_create_matrix(1, 1, &X);

  X.matrix[0][0] = 4.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_mult_number_6) {
  int res = 0;
  matrix_t A = {0};
  matrix_t Z = {0};

  s21_create_matrix(0, 4, &A);

  res = s21_mult_number(&A, 56.0, &Z);
  ck_assert_int_eq(res, ERROR);

  s21_remove_matrix(&A);
  s21_remove_matrix(&Z);
}
END_TEST

START_TEST(s21_mult_matrix_1) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(4, 4, &A);
  s21_init_matrix(1.0, &A);

  s21_create_matrix(4, 4, &B);
  s21_init_matrix(1.0, &B);

  res = s21_mult_matrix(&A, &B, &Z);
  ck_assert_int_eq(res, OK);

  s21_create_matrix(4, 4, &X);

  X.matrix[0][0] = 90.0;
  X.matrix[0][1] = 100.0;
  X.matrix[0][2] = 110.0;
  X.matrix[0][3] = 120.0;
  X.matrix[1][0] = 202.0;
  X.matrix[1][1] = 228.0;
  X.matrix[1][2] = 254.0;
  X.matrix[1][3] = 280.0;
  X.matrix[2][0] = 314.0;
  X.matrix[2][1] = 356.0;
  X.matrix[2][2] = 398.0;
  X.matrix[2][3] = 440.0;
  X.matrix[3][0] = 426.0;
  X.matrix[3][1] = 484.0;
  X.matrix[3][2] = 542.0;
  X.matrix[3][3] = 600.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_mult_matrix_2) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(4, 4, &A);
  s21_init_matrix(1.0, &A);

  s21_create_matrix(4, 4, &B);
  s21_init_matrix(-5.0, &B);

  res = s21_mult_matrix(&A, &B, &Z);
  ck_assert_int_eq(res, OK);

  s21_create_matrix(4, 4, &X);

  X.matrix[0][0] = 30.0;
  X.matrix[0][1] = 40.0;
  X.matrix[0][2] = 50.0;
  X.matrix[0][3] = 60.0;
  X.matrix[1][0] = 46.0;
  X.matrix[1][1] = 72.0;
  X.matrix[1][2] = 98.0;
  X.matrix[1][3] = 124.0;
  X.matrix[2][0] = 62.0;
  X.matrix[2][1] = 104.0;
  X.matrix[2][2] = 146.0;
  X.matrix[2][3] = 188.0;
  X.matrix[3][0] = 78.0;
  X.matrix[3][1] = 136.0;
  X.matrix[3][2] = 194.0;
  X.matrix[3][3] = 252.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_mult_matrix_3) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(4, 4, &A);
  s21_init_matrix(-5.0, &A);

  s21_create_matrix(4, 4, &B);
  s21_init_matrix(1.0, &B);

  res = s21_mult_matrix(&A, &B, &Z);
  ck_assert_int_eq(res, OK);

  s21_create_matrix(4, 4, &X);

  X.matrix[0][0] = -78.0;
  X.matrix[0][1] = -92.0;
  X.matrix[0][2] = -106.0;
  X.matrix[0][3] = -120.0;
  X.matrix[1][0] = 34.0;
  X.matrix[1][1] = 36.0;
  X.matrix[1][2] = 38.0;
  X.matrix[1][3] = 40.0;
  X.matrix[2][0] = 146.0;
  X.matrix[2][1] = 164.0;
  X.matrix[2][2] = 182.0;
  X.matrix[2][3] = 200.0;
  X.matrix[3][0] = 258.0;
  X.matrix[3][1] = 292.0;
  X.matrix[3][2] = 326.0;
  X.matrix[3][3] = 360.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_mult_matrix_4) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(4, 2, &A);
  s21_init_matrix(-5.0, &A);

  s21_create_matrix(2, 4, &B);
  s21_init_matrix(1.0, &B);

  res = s21_mult_matrix(&A, &B, &Z);
  ck_assert_int_eq(res, OK);

  s21_create_matrix(4, 4, &X);

  X.matrix[0][0] = -25.0;
  X.matrix[0][1] = -34.0;
  X.matrix[0][2] = -43.0;
  X.matrix[0][3] = -52.0;
  X.matrix[1][0] = -13.0;
  X.matrix[1][1] = -18.0;
  X.matrix[1][2] = -23.0;
  X.matrix[1][3] = -28.0;
  X.matrix[2][0] = -1.0;
  X.matrix[2][1] = -2.0;
  X.matrix[2][2] = -3.0;
  X.matrix[2][3] = -4.0;
  X.matrix[3][0] = 11.0;
  X.matrix[3][1] = 14.0;
  X.matrix[3][2] = 17.0;
  X.matrix[3][3] = 20.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_mult_matrix_5) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(1, 1, &A);
  s21_init_matrix(-5.0, &A);

  s21_create_matrix(1, 1, &B);
  s21_init_matrix(1.0, &B);

  res = s21_mult_matrix(&A, &B, &Z);
  ck_assert_int_eq(res, OK);

  s21_create_matrix(1, 1, &X);

  X.matrix[0][0] = -5.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_mult_matrix_6) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t Z = {0};

  res = s21_mult_matrix(&A, &B, &Z);
  ck_assert_int_eq(res, ERROR);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&Z);
}
END_TEST

START_TEST(s21_mult_matrix_7) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t Z = {0};

  s21_create_matrix(2, 3, &A);
  s21_init_matrix(-5.0, &A);

  s21_create_matrix(2, 2, &B);
  s21_init_matrix(1.0, &B);

  res = s21_mult_matrix(&A, &B, &Z);
  ck_assert_int_eq(res, ERROR_CALC);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&Z);
}
END_TEST

START_TEST(simple_mult) {
  const int rows = 2;
  const int cols = 3;
  matrix_t m = {0};
  s21_create_matrix(rows, cols, &m);
  matrix_t mtx = {0};
  s21_create_matrix(cols, rows, &mtx);

  for (int i = 0, c = 1; i < rows; i++)
    for (int j = 0; j < cols; j++) m.matrix[i][j] = c++;

  for (int i = 0, c = 7; i < cols; i++)
    for (int j = 0; j < rows; j++) mtx.matrix[i][j] = c++;

  matrix_t check = {0};
  s21_create_matrix(m.rows, mtx.columns, &check);
  check.matrix[0][0] = 58;
  check.matrix[0][1] = 64;
  check.matrix[1][0] = 139;
  check.matrix[1][1] = 154;

  matrix_t res = {0};
  ck_assert_int_eq(s21_mult_matrix(&m, &mtx, &res), OK);
  ck_assert_int_eq(s21_eq_matrix(&check, &res), SUCCESS);

  s21_remove_matrix(&m);
  s21_remove_matrix(&mtx);
  s21_remove_matrix(&res);
  s21_remove_matrix(&check);
}
END_TEST

START_TEST(mult_matrix2) {
  int rows = rand() % 100 + 1;
  int cols = rand() % 100 + 1;
  rows = -rows;
  cols = -cols;
  matrix_t m = {0};
  s21_create_matrix(rows, cols, &m);
  matrix_t mtx = {0};
  s21_create_matrix(cols, rows, &mtx);

  matrix_t check = {0};
  s21_create_matrix(m.rows, mtx.columns, &check);

  matrix_t res = {0};
  ck_assert_int_eq(s21_mult_matrix(&m, &mtx, &res), ERROR);

  s21_remove_matrix(&m);
  s21_remove_matrix(&mtx);
  s21_remove_matrix(&res);
  s21_remove_matrix(&check);
}
END_TEST

START_TEST(mult_matrix3) {
  matrix_t m = {0};
  int rows = 2;
  int cols = 3;
  s21_create_matrix(rows, cols, &m);
  matrix_t mtx = {0};
  int rows1 = 4;
  int cols1 = 5;
  s21_create_matrix(rows1, cols1, &mtx);

  matrix_t check = {0};
  s21_create_matrix(m.rows, mtx.columns, &check);

  matrix_t res = {0};
  ck_assert_int_eq(s21_mult_matrix(&m, &mtx, &res), ERROR_CALC);

  s21_remove_matrix(&m);
  s21_remove_matrix(&mtx);
  s21_remove_matrix(&res);
  s21_remove_matrix(&check);
}
END_TEST

START_TEST(s21_transpose_1) {
  int res = 0;
  matrix_t A = {0};

  res = s21_transpose(&A, &A);
  ck_assert_int_eq(res, ERROR);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_transpose_2) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};

  s21_create_matrix(3, 4, &A);
  s21_init_matrix(1.0, &A);
  res = s21_transpose(&A, &B);
  ck_assert_int_eq(res, OK);
  res = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(res, FAILURE);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(s21_transpose_3) {
  int res = 0;
  matrix_t A = {0};
  matrix_t T = {0};
  matrix_t B = {0};

  s21_create_matrix(2, 3, &A);
  s21_init_matrix(1.0, &A);
  s21_create_matrix(3, 2, &B);
  B.matrix[0][0] = 1.0;
  B.matrix[0][1] = 4.0;
  B.matrix[1][0] = 2.0;
  B.matrix[1][1] = 5.0;
  B.matrix[2][0] = 3.0;
  B.matrix[2][1] = 6.0;

  res = s21_transpose(&A, &T);
  ck_assert_int_eq(res, OK);
  res = s21_eq_matrix(&T, &B);
  ck_assert_int_eq(res, SUCCESS);
  s21_remove_matrix(&A);
  s21_remove_matrix(&T);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(transpose_matrix2) {
  int rows = rand() % 100 + 1;
  int cols = rand() % 100 + 1;
  rows = -rows;
  cols = -cols;
  matrix_t m = {0};
  s21_create_matrix(rows, cols, &m);
  matrix_t check = {0};
  s21_create_matrix(cols, rows, &check);

  matrix_t res = {0};
  ck_assert_int_eq(s21_transpose(&m, &res), ERROR);

  s21_remove_matrix(&m);
  s21_remove_matrix(&res);
  s21_remove_matrix(&check);
}
END_TEST

START_TEST(s21_determinant_1) {
  int res = 0;
  double determ = 0.0;
  matrix_t A = {NULL, 0, 0};

  res = s21_determinant(&A, &determ);
  ck_assert_int_eq(res, ERROR);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_determinant_2) {
  int res = 0;
  double determ = 0.0;
  matrix_t A = {0};
  s21_create_matrix(1, 2, &A);

  res = s21_determinant(&A, &determ);
  ck_assert_int_eq(res, ERROR_CALC);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_determinant_3) {
  double determ = 0.0;
  matrix_t A = {0};

  s21_create_matrix(1, 1, &A);
  s21_init_matrix(1.0, &A);

  s21_determinant(&A, &determ);
  ck_assert_double_eq(determ, 1.0);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_determinant_4) {
  double determ = 0.0;
  matrix_t A = {0};

  s21_create_matrix(3, 3, &A);
  A.matrix[0][0] = 4000.87;
  A.matrix[0][1] = 454.0;
  A.matrix[0][2] = 4.0;
  A.matrix[1][0] = 2.0;
  A.matrix[1][1] = 1.897;
  A.matrix[1][2] = 3.0;
  A.matrix[2][0] = 0.000087;
  A.matrix[2][1] = 2.5668;
  A.matrix[2][2] = 1.0;

  s21_determinant(&A, &determ);
  ck_assert_double_eq(determ, -24105.996724156);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_determinant_5) {
  double determ = 0.0;
  matrix_t A = {0};

  s21_create_matrix(2, 2, &A);
  A.matrix[0][0] = -677700.0;
  A.matrix[0][1] = 654.0;
  A.matrix[1][0] = 2.0;
  A.matrix[1][1] = -0.00001;

  s21_determinant(&A, &determ);
  ck_assert_double_eq(determ, -1301.223);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_determinant_6) {
  double determ = 0.0;
  matrix_t A = {0};

  s21_create_matrix(8, 8, &A);
  s21_init_matrix(1.0, &A);

  s21_determinant(&A, &determ);
  ck_assert_double_eq(determ, 0.0);
  s21_remove_matrix(&A);
}
END_TEST
START_TEST(determinant1) {
  const int size = 5;
  matrix_t m = {0};
  s21_create_matrix(size, size, &m);

  for (int i = 0; i < size; i++)
    for (int j = 0; j < size; j++) m.matrix[i][j] = j;
  double res = 0;
  int code = s21_determinant(&m, &res);
  ck_assert_double_eq_tol(res, 0, 1e-6);
  ck_assert_int_eq(code, OK);

  s21_remove_matrix(&m);
}
END_TEST

START_TEST(determinant2) {
  const int size = 4;
  matrix_t m = {0};
  s21_create_matrix(size, size, &m);

  for (int i = 0; i < size; i++)
    for (int j = 0; j < size; j++) m.matrix[i][j] = j + i;

  double res = 0;
  int code = s21_determinant(&m, &res);
  ck_assert_double_eq_tol(res, 0, 1e-6);
  ck_assert_int_eq(code, OK);

  s21_remove_matrix(&m);
}
END_TEST

START_TEST(determinant3) {
  const int size = 5;
  matrix_t m = {0};
  s21_create_matrix(size, size, &m);
  m.matrix[0][1] = 6;
  m.matrix[0][2] = -2;
  m.matrix[0][3] = -1;
  m.matrix[0][4] = 5;
  m.matrix[1][3] = -9;
  m.matrix[1][4] = -7;
  m.matrix[2][1] = 15;
  m.matrix[2][2] = 35;
  m.matrix[3][1] = -1;
  m.matrix[3][2] = -11;
  m.matrix[3][3] = -2;
  m.matrix[3][4] = 1;
  m.matrix[4][0] = -2;
  m.matrix[4][1] = -2;
  m.matrix[4][2] = 3;
  m.matrix[4][4] = -2;

  double res = 0;
  int code = s21_determinant(&m, &res);
  ck_assert_double_eq_tol(res, 2480, 1e-6);
  ck_assert_int_eq(code, OK);

  s21_remove_matrix(&m);
}
END_TEST

START_TEST(determinant4) {
  const int size = 3;
  matrix_t m = {0};
  s21_create_matrix(size, size, &m);
  m.matrix[0][0] = 2;
  m.matrix[0][1] = 3;
  m.matrix[0][2] = 1;
  m.matrix[1][0] = 7;
  m.matrix[1][1] = 4;
  m.matrix[1][2] = 1;
  m.matrix[2][0] = 9;
  m.matrix[2][1] = -2;
  m.matrix[2][2] = 1;

  double res = 0;
  int code = s21_determinant(&m, &res);
  ck_assert_double_eq_tol(res, -32, 1e-6);
  ck_assert_int_eq(code, OK);

  s21_remove_matrix(&m);
}
END_TEST

START_TEST(determinant5) {
  const int size = 2;
  matrix_t m = {0};
  s21_create_matrix(size, size, &m);
  m.matrix[0][0] = -5;
  m.matrix[0][1] = -4;
  m.matrix[1][0] = -2;
  m.matrix[1][1] = -3;

  double res = 0;
  int code = s21_determinant(&m, &res);
  ck_assert_double_eq_tol(res, 7, 1e-6);
  ck_assert_int_eq(code, OK);

  s21_remove_matrix(&m);
}
END_TEST

START_TEST(determinant6) {
  const int size = 1;
  matrix_t m = {0};
  s21_create_matrix(size, size, &m);
  m.matrix[0][0] = -5;

  double res = 0;
  int code = s21_determinant(&m, &res);
  ck_assert_double_eq_tol(res, -5, 1e-6);
  ck_assert_int_eq(code, OK);

  s21_remove_matrix(&m);
}
END_TEST

START_TEST(determinant7) {
  matrix_t m = {0};
  int rows = rand() % 100 + 1;
  rows = -rows;
  s21_create_matrix(rows, rows, &m);
  double res = 0;
  int code = s21_determinant(&m, &res);
  ck_assert_int_eq(code, ERROR);

  s21_remove_matrix(&m);
}
END_TEST

START_TEST(s21_calc_complements_1) {
  int res = 0;
  matrix_t A = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(4, 4, &A);
  A.matrix[0][0] = 112.0;
  A.matrix[0][1] = 243.0;
  A.matrix[0][2] = 3.0;
  A.matrix[0][3] = -8.0;
  A.matrix[1][0] = 47.0;
  A.matrix[1][1] = 51.0;
  A.matrix[1][2] = -66.0;
  A.matrix[1][3] = 99.0;
  A.matrix[2][0] = -74.0;
  A.matrix[2][1] = 85.0;
  A.matrix[2][2] = 97.0;
  A.matrix[2][3] = 63.0;
  A.matrix[3][0] = -13.0;
  A.matrix[3][1] = 79.0;
  A.matrix[3][2] = -99.0;
  A.matrix[3][3] = -121.0;

  res = s21_calc_complements(&A, &Z);
  ck_assert_int_eq(res, OK);

  s21_create_matrix(4, 4, &X);

  X.matrix[0][0] = -2879514.0;
  X.matrix[0][1] = -1236631.0;
  X.matrix[0][2] = -1685096.0;
  X.matrix[0][3] = 880697.0;
  X.matrix[1][0] = 1162090.0;
  X.matrix[1][1] = -714015.0;
  X.matrix[1][2] = 4046255.0;
  X.matrix[1][3] = -3901600.0;
  X.matrix[2][0] = 4362897.0;
  X.matrix[2][1] = -2049432.0;
  X.matrix[2][2] = -532912.0;
  X.matrix[2][3] = -1370781.0;
  X.matrix[3][0] = 3412773.0;
  X.matrix[3][1] = -1569493.0;
  X.matrix[3][2] = 3144517.0;
  X.matrix[3][3] = 1284666.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_calc_complements_2) {
  int res = 0;
  matrix_t A = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(4, 4, &A);
  A.matrix[0][0] = -1.0;
  A.matrix[0][1] = 2.0;
  A.matrix[0][2] = 7.0;
  A.matrix[0][3] = 9.0;
  A.matrix[1][0] = 1.0;
  A.matrix[1][1] = 0.0;
  A.matrix[1][2] = 0.0;
  A.matrix[1][3] = 0.0;
  A.matrix[2][0] = 47.0;
  A.matrix[2][1] = 13.0;
  A.matrix[2][2] = 17.0;
  A.matrix[2][3] = 21.0;
  A.matrix[3][0] = 22.0;
  A.matrix[3][1] = 7.0;
  A.matrix[3][2] = 1.0;
  A.matrix[3][3] = 3.0;

  res = s21_calc_complements(&A, &Z);
  ck_assert_int_eq(res, OK);

  s21_create_matrix(4, 4, &X);

  X.matrix[0][0] = 0.0;
  X.matrix[0][1] = -30.0;
  X.matrix[0][2] = -108.0;
  X.matrix[0][3] = 106.0;
  X.matrix[1][0] = 138.0;
  X.matrix[1][1] = -726.0;
  X.matrix[1][2] = -1137.0;
  X.matrix[1][3] = 1061.0;
  X.matrix[2][0] = 0.0;
  X.matrix[2][1] = 12.0;
  X.matrix[2][2] = 57.0;
  X.matrix[2][3] = -47.0;
  X.matrix[3][0] = 0.0;
  X.matrix[3][1] = 6.0;
  X.matrix[3][2] = -75.0;
  X.matrix[3][3] = 57.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_calc_complements_3) {
  int res = 0;
  matrix_t A = {0};
  matrix_t Z = {0};

  s21_create_matrix(4, 3, &A);
  A.matrix[0][0] = -1.0;
  A.matrix[0][1] = 2.0;
  A.matrix[0][2] = 7.0;
  A.matrix[1][0] = 1.0;
  A.matrix[1][1] = 0.0;
  A.matrix[1][2] = 0.0;
  A.matrix[2][0] = 47.0;
  A.matrix[2][1] = 13.0;
  A.matrix[2][2] = 17.0;
  A.matrix[3][0] = 22.0;
  A.matrix[3][1] = 7.0;
  A.matrix[3][2] = 1.0;

  res = s21_calc_complements(&A, &Z);

  s21_remove_matrix(&A);

  ck_assert_int_eq(res, ERROR_CALC);
}
END_TEST

START_TEST(s21_calc_complements_4) {
  int res = 0;
  matrix_t A = {0};
  matrix_t Z = {0};

  s21_create_matrix(2, 2, &A);
  s21_remove_matrix(&A);  // NULL

  res = s21_calc_complements(&A, &Z);
  ck_assert_int_eq(res, ERROR);
}
END_TEST

START_TEST(s21_inverse_matrix_1) {
  int res = 0;
  matrix_t A = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(4, 4, &A);
  A.matrix[0][0] = -1.0;
  A.matrix[0][1] = 2.0;
  A.matrix[0][2] = 7.0;
  A.matrix[0][3] = 9.0;
  A.matrix[1][0] = 1.0;
  A.matrix[1][1] = 0.0;
  A.matrix[1][2] = 0.0;
  A.matrix[1][3] = 0.0;
  A.matrix[2][0] = 47.0;
  A.matrix[2][1] = 13.0;
  A.matrix[2][2] = 17.0;
  A.matrix[2][3] = 21.0;
  A.matrix[3][0] = 22.0;
  A.matrix[3][1] = 7.0;
  A.matrix[3][2] = 1.0;
  A.matrix[3][3] = 3.0;

  s21_inverse_matrix(&A, &Z);

  s21_create_matrix(4, 4, &X);

  X.matrix[0][0] = 0.0;
  X.matrix[0][1] = 1.0;
  X.matrix[0][2] = 0.0;
  X.matrix[0][3] = 0.0;
  X.matrix[1][0] = -5.0 / 23.0;
  X.matrix[1][1] = -121.0 / 23.0;
  X.matrix[1][2] = 2.0 / 23.0;
  X.matrix[1][3] = 1.0 / 23.0;
  X.matrix[2][0] = -18.0 / 23.0;
  X.matrix[2][1] = -379.0 / 46.0;
  X.matrix[2][2] = 19.0 / 46.0;
  X.matrix[2][3] = -25.0 / 46.0;
  X.matrix[3][0] = 53.0 / 69.0;
  X.matrix[3][1] = 1061.0 / 138.0;
  X.matrix[3][2] = -47.0 / 138.0;
  X.matrix[3][3] = 19.0 / 46.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_inverse_matrix_2) {
  int res;
  matrix_t A = {0};
  matrix_t Z = {0};

  s21_create_matrix(1, 1, &A);
  A.matrix[0][0] = 2;

  res = s21_inverse_matrix(&A, &Z);

  ck_assert_double_eq(1.0 / A.matrix[0][0], 0.5);
  ck_assert_int_eq(res, OK);

  s21_remove_matrix(&A);
  s21_remove_matrix(&Z);
}
END_TEST

START_TEST(s21_inverse_matrix_3) {
  int res = 0;
  matrix_t A = {0};
  matrix_t Z = {0};

  s21_create_matrix(3, 3, &A);
  s21_init_matrix(1.0, &A);

  res = s21_inverse_matrix(&A, &Z);
  ck_assert_int_eq(res, ERROR_CALC);
  s21_remove_matrix(&Z);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_inverse_matrix_4) {
  int res = 0;
  matrix_t A = {0};
  matrix_t Z = {0};

  s21_create_matrix(3, 2, &A);
  s21_init_matrix(1.0, &A);

  res = s21_inverse_matrix(&A, &Z);
  ck_assert_int_eq(res, ERROR_CALC);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_inverse_matrix_5) {
  int res = 0;
  matrix_t A = {0};
  matrix_t Z = {0};

  s21_create_matrix(3, 3, &A);
  s21_init_matrix(1.0, &A);
  s21_remove_matrix(&A);

  res = s21_inverse_matrix(&A, &Z);

  ck_assert_int_eq(res, ERROR);
}
END_TEST

START_TEST(s21_inverse_matrix_6) {
  int res = 0;
  matrix_t A = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(3, 3, &A);
  A.matrix[0][0] = 1.0;
  A.matrix[0][1] = 2.0;
  A.matrix[0][2] = 3.0;
  A.matrix[1][0] = 4.0;
  A.matrix[1][1] = 5.0;
  A.matrix[1][2] = 6.0;
  A.matrix[2][0] = 8.0;
  A.matrix[2][1] = 6.0;
  A.matrix[2][2] = 9.0;

  s21_inverse_matrix(&A, &Z);

  s21_create_matrix(3, 3, &X);

  X.matrix[0][0] = -3.0 / 5.0;
  X.matrix[0][1] = 0.0;
  X.matrix[0][2] = 1.0 / 5.0;
  ;
  X.matrix[1][0] = -4.0 / 5.0;
  X.matrix[1][1] = 1.0;
  X.matrix[1][2] = -2.0 / 5.0;
  X.matrix[2][0] = 16.0 / 15.0;
  X.matrix[2][1] = -2.0 / 3.0;
  X.matrix[2][2] = 1.0 / 5.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_inverse_1) {
  matrix_t A, C;
  s21_create_matrix(3, 3, &A);
  s21_create_matrix(3, 3, &C);
  C.matrix[0][0] = 1.0;
  C.matrix[0][1] = -1.0;
  C.matrix[0][2] = 1.0;
  C.matrix[1][0] = -38.0;
  C.matrix[1][1] = 41.0;
  C.matrix[1][2] = -34.0;
  C.matrix[2][0] = 27.0;
  C.matrix[2][1] = -29.0;
  C.matrix[2][2] = 24.0;
  A.matrix[0][0] = 2.0;
  A.matrix[0][1] = 5.0;
  A.matrix[0][2] = 7.0;
  A.matrix[1][0] = 6.0;
  A.matrix[1][1] = 3.0;
  A.matrix[1][2] = 4.0;
  A.matrix[2][0] = 5.0;
  A.matrix[2][1] = -2.0;
  A.matrix[2][2] = -3.0;
  matrix_t B;
  s21_inverse_matrix(&A, &B);
  int res = s21_eq_matrix(&B, &C);
  ck_assert_int_eq(res, 1);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&C);
}
END_TEST

START_TEST(test_one_by_one) {
  matrix_t m = {0};
  matrix_t result = {0};
  int codec = s21_create_matrix(1, 1, &m);
  if (!codec) {
    m.matrix[0][0] = 1431.12312331;
    int code = s21_inverse_matrix(&m, &result);
    ck_assert_int_eq(1.0 / m.matrix[0][0], 0);
    ck_assert_int_eq(code, OK);
    s21_remove_matrix(&m);
    s21_remove_matrix(&result);
  }
}
END_TEST

START_TEST(test_zero_det) {
  matrix_t m = {0};
  matrix_t result = {0};
  int codec = s21_create_matrix(1, 1, &m);
  if (!codec) {
    int code = s21_inverse_matrix(&m, &result);
    ck_assert_int_eq(code, ERROR_CALC);
    s21_remove_matrix(&m);
    s21_remove_matrix(&result);
  }
}
END_TEST

START_TEST(test_incorrect) {
  matrix_t m = {0};
  matrix_t result = {0};
  int code = s21_inverse_matrix(&m, &result);
  ck_assert_int_eq(code, ERROR);
}
END_TEST

START_TEST(determinant) {
  const int size = 2;
  matrix_t m = {0};
  s21_create_matrix(size, size, &m);
  m.matrix[0][0] = 1;
  m.matrix[0][1] = 1;
  m.matrix[1][0] = 1;
  m.matrix[1][1] = 1;

  matrix_t result = {0};
  int code = s21_inverse_matrix(&m, &result);
  ck_assert_int_eq(code, ERROR_CALC);

  s21_remove_matrix(&m);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_not_sqare) {
  matrix_t m = {0};
  matrix_t result = {0};
  int codec = s21_create_matrix(1, 4, &m);
  if (!codec) {
    int code = s21_inverse_matrix(&m, &result);
    ck_assert_int_eq(code, ERROR_CALC);
    s21_remove_matrix(&m);
  }
}
END_TEST

START_TEST(test_normal) {
  matrix_t m = {0};
  matrix_t expected = {0};
  int codec1, codec2;
  codec1 = s21_create_matrix(3, 3, &m);
  if (!codec1) codec2 = s21_create_matrix(3, 3, &expected);

  if (!codec1 && !codec2) {
    m.matrix[0][0] = 2;
    m.matrix[0][1] = 5;
    m.matrix[0][2] = 7;

    m.matrix[1][0] = 6;
    m.matrix[1][1] = 3;
    m.matrix[1][2] = 4;

    m.matrix[2][0] = 5;
    m.matrix[2][1] = -2;
    m.matrix[2][2] = -3;

    expected.matrix[0][0] = 1;
    expected.matrix[0][1] = -1;
    expected.matrix[0][2] = 1;

    expected.matrix[1][0] = -38;
    expected.matrix[1][1] = 41;
    expected.matrix[1][2] = -34;

    expected.matrix[2][0] = 27;
    expected.matrix[2][1] = -29;
    expected.matrix[2][2] = 24;
    matrix_t result = {0};
    int code = s21_inverse_matrix(&m, &result);

    ck_assert_int_eq(s21_eq_matrix(&result, &expected), SUCCESS);
    ck_assert_int_eq(code, OK);

    s21_remove_matrix(&m);
    s21_remove_matrix(&result);
    s21_remove_matrix(&expected);
  }
}
END_TEST

START_TEST(inverse) {
  const int size = 3;
  matrix_t m = {0};
  s21_create_matrix(size, size, &m);

  m.matrix[0][0] = 2;
  m.matrix[0][1] = 5;
  m.matrix[0][2] = 7;
  m.matrix[1][0] = 6;
  m.matrix[1][1] = 3;
  m.matrix[1][2] = 4;
  m.matrix[2][0] = 5;
  m.matrix[2][1] = -2;
  m.matrix[2][2] = -3;

  matrix_t res = {0};
  s21_inverse_matrix(&m, &res);

  matrix_t expected = {0};
  s21_create_matrix(size, size, &expected);
  expected.matrix[0][0] = 1;
  expected.matrix[0][1] = -1;
  expected.matrix[0][2] = 1;
  expected.matrix[1][0] = -38;
  expected.matrix[1][1] = 41;
  expected.matrix[1][2] = -34;
  expected.matrix[2][0] = 27;
  expected.matrix[2][1] = -29;
  expected.matrix[2][2] = 24;

  ck_assert_int_eq(s21_eq_matrix(&expected, &res), SUCCESS);

  s21_remove_matrix(&expected);
  s21_remove_matrix(&res);
  s21_remove_matrix(&m);
}
END_TEST

Suite *s21_create_matrix_suite(void) {
  Suite *s;
  TCase *tc;

  s = suite_create("\033[42ms21_create_matrix\033[0m\033[0m");

  tc = tcase_create("\033[42mcase_create_matrix\033[0m\033[0m");

  tcase_add_test(tc, s21_create_matrix_1);
  tcase_add_test(tc, s21_create_matrix_2);
  tcase_add_test(tc, s21_create_matrix_3);
  tcase_add_test(tc, s21_create_matrix_4);
  tcase_add_test(tc, s21_create_matrix_5);
  tcase_add_test(tc, create_1);
  tcase_add_test(tc, create_2);
  tcase_add_test(tc, create_3);

  suite_add_tcase(s, tc);
  return s;
}

Suite *s21_eq_matrix_suite(void) {
  Suite *s;
  TCase *tc;

  s = suite_create("\033[42ms21_eq_matrix\033[0m");

  tc = tcase_create("\033[42mcase_eq_matrix\033[0m");

  tcase_add_test(tc, s21_eq_matrix_1);
  tcase_add_test(tc, s21_eq_matrix_2);
  tcase_add_test(tc, s21_eq_matrix_3);
  tcase_add_test(tc, s21_eq_matrix_4);
  tcase_add_test(tc, s21_eq_matrix_5);
  tcase_add_test(tc, zero_matrix);
  tcase_add_test(tc, zero_matrix_1);
  tcase_add_test(tc, casual_matrix_1);
  tcase_add_test(tc, casual_matrix_2);
  tcase_add_test(tc, casual_matrix_3);
  tcase_add_test(tc, casual_matrix_4);
  tcase_add_test(tc, casual_matrix_5);
  tcase_add_test(tc, casual_matrix_6);
  tcase_add_test(tc, casual_matrix_7);
  tcase_add_test(tc, casual_matrix_8);

  suite_add_tcase(s, tc);
  return s;
}

Suite *s21_sum_matrix_suite(void) {
  Suite *s;
  TCase *tc;

  s = suite_create("\033[42ms21_sum_matrix\033[0m");

  tc = tcase_create("\033[42mcase_sum_matrix\033[0m");

  tcase_add_test(tc, s21_sum_matrix_1);
  tcase_add_test(tc, s21_sum_matrix_2);
  tcase_add_test(tc, s21_sum_matrix_3);
  tcase_add_test(tc, s21_sum_matrix_4);
  tcase_add_loop_test(tc, sum_matrix2, 0, 10);

  suite_add_tcase(s, tc);
  return s;
}

Suite *s21_sub_matrix_suite(void) {
  Suite *s;
  TCase *tc;

  s = suite_create("\033[42ms21_sub_matrix\033[0m");

  tc = tcase_create("\033[42mcase_sub_matrix\033[0m");

  tcase_add_test(tc, s21_sub_matrix_1);
  tcase_add_test(tc, s21_sub_matrix_2);
  tcase_add_test(tc, s21_sub_matrix_3);
  tcase_add_test(tc, s21_sub_matrix_4);
  tcase_add_loop_test(tc, sub_matrix2, 0, 10);

  suite_add_tcase(s, tc);
  return s;
}

Suite *s21_mult_number_matrix_suite(void) {
  Suite *s;
  TCase *tc;

  s = suite_create("\033[42ms21_mult_number_matrix\033[0m");

  tc = tcase_create("\033[42mcase_mult_number_matrix\033[0m");

  tcase_add_test(tc, s21_mult_number_1);
  tcase_add_test(tc, s21_mult_number_2);
  tcase_add_test(tc, s21_mult_number_3);
  tcase_add_test(tc, s21_mult_number_4);
  tcase_add_test(tc, s21_mult_number_5);
  tcase_add_test(tc, s21_mult_number_6);
  tcase_add_loop_test(tc, mult_matrix2, 0, 10);
  tcase_add_test(tc, mult_matrix3);
  tcase_add_test(tc, simple_mult);

  suite_add_tcase(s, tc);
  return s;
}

Suite *s21_mult_matrix_suite(void) {
  Suite *s;
  TCase *tc;

  s = suite_create("\033[42ms21_mult_matrix\033[0m");

  tc = tcase_create("\033[42mcase_mult_matrix\033[0m");

  tcase_add_test(tc, s21_mult_matrix_1);
  tcase_add_test(tc, s21_mult_matrix_2);
  tcase_add_test(tc, s21_mult_matrix_3);
  tcase_add_test(tc, s21_mult_matrix_4);
  tcase_add_test(tc, s21_mult_matrix_5);
  tcase_add_test(tc, s21_mult_matrix_6);
  tcase_add_test(tc, s21_mult_matrix_7);

  suite_add_tcase(s, tc);
  return s;
}

Suite *s21_transpose_suite(void) {
  Suite *s;
  TCase *tc;

  s = suite_create("\033[42ms21_transpose\033[0m");

  tc = tcase_create("\033[42mcase_transpose\033[0m");

  tcase_add_test(tc, s21_transpose_1);
  tcase_add_test(tc, s21_transpose_2);
  tcase_add_test(tc, s21_transpose_3);
  tcase_add_loop_test(tc, transpose_matrix2, 0, 10);

  suite_add_tcase(s, tc);
  return s;
}

Suite *s21_determinant_suite(void) {
  Suite *s;
  TCase *tc;

  s = suite_create("\033[42ms21_determinant\033[0m");

  tc = tcase_create("\033[42mcase_determinant\033[0m");

  tcase_add_test(tc, s21_determinant_1);
  tcase_add_test(tc, s21_determinant_2);
  tcase_add_test(tc, s21_determinant_3);
  tcase_add_test(tc, s21_determinant_4);
  tcase_add_test(tc, s21_determinant_5);
  tcase_add_test(tc, s21_determinant_6);
  tcase_add_test(tc, determinant1);
  tcase_add_test(tc, determinant2);
  tcase_add_test(tc, determinant3);
  tcase_add_test(tc, determinant4);
  tcase_add_test(tc, determinant5);
  tcase_add_test(tc, determinant6);
  tcase_add_loop_test(tc, determinant7, 0, 10);

  suite_add_tcase(s, tc);
  return s;
}

Suite *s21_calc_complements_suite(void) {
  Suite *s;
  TCase *tc;

  s = suite_create("\033[42ms21_calc_complements\033[0m");

  tc = tcase_create("\033[42mcase_calc_complements\033[0m");

  tcase_add_test(tc, s21_calc_complements_1);
  tcase_add_test(tc, s21_calc_complements_2);
  tcase_add_test(tc, s21_calc_complements_3);
  tcase_add_test(tc, s21_calc_complements_4);
  tcase_add_test(tc, test_normal);
  tcase_add_test(tc, test_not_sqare);
  tcase_add_test(tc, test_incorrect);
  tcase_add_test(tc, test_one_by_one);

  suite_add_tcase(s, tc);
  return s;
}

Suite *s21_inverse_matrix_suite(void) {
  Suite *s;
  TCase *tc;

  s = suite_create("\033[42ms21_inverse_matrix\033[0m");

  tc = tcase_create("\033[42mcase_inverse_matrix\033[0m");

  tcase_add_test(tc, s21_inverse_matrix_1);
  tcase_add_test(tc, s21_inverse_matrix_2);
  tcase_add_test(tc, s21_inverse_matrix_3);
  tcase_add_test(tc, s21_inverse_matrix_4);
  tcase_add_test(tc, s21_inverse_matrix_5);
  tcase_add_test(tc, s21_inverse_matrix_6);
  tcase_add_test(tc, inverse);
  tcase_add_test(tc, s21_inverse_1);
  tcase_add_test(tc, determinant);
  tcase_add_test(tc, test_normal);
  tcase_add_test(tc, test_not_sqare);
  tcase_add_test(tc, test_zero_det);
  tcase_add_test(tc, test_incorrect);
  tcase_add_test(tc, test_one_by_one);

  suite_add_tcase(s, tc);
  return s;
}

void case_test(Suite *s, int *fail) {
  SRunner *runner = srunner_create(s);
  srunner_set_fork_status(runner, CK_NOFORK);
  srunner_run_all(runner, CK_NORMAL);
  *fail = srunner_ntests_failed(runner);
  srunner_free(runner);
}

int main() {
  int fail = 0;

  case_test(s21_create_matrix_suite(), &fail);
  case_test(s21_eq_matrix_suite(), &fail);
  case_test(s21_sum_matrix_suite(), &fail);
  case_test(s21_sub_matrix_suite(), &fail);
  case_test(s21_mult_number_matrix_suite(), &fail);
  case_test(s21_mult_matrix_suite(), &fail);
  case_test(s21_transpose_suite(), &fail);
  case_test(s21_determinant_suite(), &fail);
  case_test(s21_calc_complements_suite(), &fail);
  case_test(s21_inverse_matrix_suite(), &fail);

  return fail;
}
