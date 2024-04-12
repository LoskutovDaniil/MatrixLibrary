#ifndef S21_MATRIX_H
#define S21_MATRIX_H

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define OK 0  // для всего, кроме сравнений матриц
#define ERROR 1  // для всего, кроме сравнений матриц
#define ERROR_CALC 2  // для всего, кроме сравнений матриц

#define SUCCESS 1  // для сравнений матриц
#define FAILURE 0  // для сравнений матриц

typedef struct matrix_struct {
  double **matrix;

  int rows;     // строки
  int columns;  // столбцы

} matrix_t;

// создание матриц
int s21_create_matrix(int rows, int columns, matrix_t *result);
// очистка матриц
void s21_remove_matrix(matrix_t *A);
// сравнение матриц
int s21_eq_matrix(matrix_t *A, matrix_t *B);
// сложение матриц
int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
// вычитание матриц
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
// Умножение матрицы на число
int s21_mult_number(matrix_t *A, double number, matrix_t *result);
// Умножение двух матриц
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
// Транспонирование матрицы
int s21_transpose(matrix_t *A, matrix_t *result);
// Минор матрицы и матрица алгебраических дополнений
int s21_calc_complements(matrix_t *A, matrix_t *result);
// Определитель матрицы
int s21_determinant(matrix_t *A, double *result);
// Обратная матрица
int s21_inverse_matrix(matrix_t *A, matrix_t *result);

void support_print_matrix(matrix_t *result);
void support_free_matrix(matrix_t *result);
int support_eq_size_matrix(matrix_t *A, matrix_t *B);
void support_search_minor(matrix_t *A, matrix_t *minor, int pass_column,
                          int pass_row);
int support_be_matrix(matrix_t *A);

#endif